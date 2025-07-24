%% gozaresh kar :
% get fft of the signal
%

clc;clear all; close all;
% DBIM + FDTD Simulation of Anthropomorphic Two-layer Calcaneus Phantom

%% 1. Domain and Parameters

f0 = 1e9;
c=3e8;
landa=c/f0;
w=2*pi*f0;

Nx = 200; Ny = 200;      % grid points
% dx = 1e-3; dy = dx;
dx = 0.2*landa/(10*sqrt(30)); dy = dx;   % spatial step (1 mm)   dx<landa/(10 sqrt(er_max))   dx<
Nt = 4000;               % time steps
start_n= 1500;
c = 3e8;
dt = dx/(sqrt(2*8)*c);         % CFL stability


eps0 = 8.85e-12;
mu0 = 4*pi*1e-7;

Nt_sample=Nt/10;


% source
 t0 =1e-9;      % Gaussian pulse center and delay
t = (0:Nt-1)*dt;
 % source_pulse = 10*exp(-((t-t0)/(t0/4)).^2).*cos(2*pi*f0*t);
 source_pulse = 2*sin(w*t) ;



%pml
pml_L=60;
sigma_max=20;
R_out=20 % edit kon function ra
R_antenna=32;
holdOffset=(Nx/2-pml_L)-R_out;
offset=holdOffset-2;
nx_new=Nx-2*(pml_L+offset);
ny_new=Ny-2*(pml_L+offset);
% nx_new=22;
% ny_new=22;

Limit_x=pml_L+offset;

%% 2. Phantom Definition

[ eps_r_bg,er_inf_cortical, er_inf_normalTrab,eps_r_true] = define_phantom(Nx, Ny,pml_L,w);

documentation.define_phantom = @define_phantom;

%% 3. Antenna Array Configuration
Ntx = 16;               % transmissions (or 24)
[tx_idx, ty_idx] = circular_array(Nx, Ny, Ntx,pml_L);

E_inc=zeros(Ntx,nx_new,ny_new);
Sij_true=zeros(Ntx,Ntx);
Sij_inc=zeros(Ntx,Ntx);
documentation.circular_array = @circular_array;

%% 4. Precompute Incident Fields (Background)   Edit this


er_test=ones(Nx,Ny)*real(eps_r_bg);
sigma_bg=abs(imag(eps_r_bg)*w)*eps0;           sigma_test=ones(Nx,Ny)*sigma_bg;


% Ek_inc = zeros(k,Nx,Ny);
% Sij_inc= zeros (k,k)
for k = 1:Ntx
    [E_inc(k,:,:) , Sij_inc(k,:)] = fdtd_simulation(Limit_x,nx_new,ny_new,Nx, Ny, dx, dy, Nt, er_test, sigma_test,k, tx_idx, ty_idx, source_pulse,pml_L,sigma_max,sigma_bg,dt,w);
end

% S_true
er_original= real(eps_r_true);
sigma_original= abs(w*imag(eps_r_true))*eps0;

for k = 1:Ntx
    [ trash , Sij_true(k,:)] = fdtd_simulation(Limit_x,nx_new,ny_new,Nx, Ny, dx, dy, Nt, er_original, sigma_original,k, tx_idx, ty_idx, source_pulse,pml_L,sigma_max,sigma_bg,dt,w);
end
clear trash;




%% born

%% 5. DBIM + IMATCS Inversion Loop
%next  : formul 13-2 paynname
% Sij_predected =Sij_inc+ ?
Holder_Sij=0;
ds=dx*dy;
Holder_Integrand=0;

eps_r_predect=real(eps_r_bg) * ones(nx_new,ny_new);

ai=1;
aj=1;
% ============================== calculate the Sij scatered : =======================================

% build the Sensivity matrix A-----------------

eps_r_new =er_test;
sigma_new =sigma_test;
% eps_r_new( Limit_x + 1 : Nx-Limit_x,Limit_x+1 :Nx-Limit_x )= eps_r_new(Limit_x+1:Nx-Limit_x ,Limit_x+1 :Nx-Limit_x) +1  ;

Np = nx_new*ny_new;          % total number of pixels
A= zeros(Ntx^2, Np);  % 256 x Np
for itr=1:40 % iteration

    %first guess
    if itr>=1
        % forward solver - FDTD
        % er_test=ones(Nx,Ny)*eps_r_bg;
        % sigma_test=ones(Nx,Ny)*sigma_bg;
        
        for k = 1:Ntx
            [Einc_k(k,:,:) , Sij_inc_k(k,:)] = fdtd_simulation(Limit_x,nx_new,ny_new,Nx, Ny, dx, dy, Nt, eps_r_new, sigma_new,k, tx_idx, ty_idx, source_pulse,pml_L,sigma_max,sigma_bg,dt,w);
        end

    else
        % at first guss
        Einc_k=E_inc;
        Sij_inc_k=Sij_inc;
    end
    A=Calc_Jacobi(Np,Ntx,Einc_k,ds);

    %edit this section !!
    Sij_predected=zeros(Ntx,Ntx);
    %first guess :
    for(I=1:Ntx)
        for(J=1:Ntx)

            Holder_Integrand=0;
            for Ix=1:nx_new
                for Iy=1:ny_new
                    contrast = eps_r_predect(Ix,Iy)-eps_r_bg;
                    hold=Einc_k(I,Ix,Iy).*Einc_k(J,Ix,Iy);
                    A_size=size(hold);
                    Holder_Integrand=Holder_Integrand+hold*ds *contrast;
                end
            end
            Sij_predected(I,J) = Sij_inc_k(I,J)-1i*(w*eps0)/(2*ai*aj).*Holder_Integrand;
        end
    end

    % compute the error
    Sij_error=Sij_true-Sij_predected;
    sij_holder(itr,:,:)=Sij_predected; % for comparing if we have any change ...

    DeltaS_vec = reshape(Sij_error.', [], 1);  % column vector, size 256x1

    %4. Solve the Inverse Problem (Regularized)
    lambda = 1e-5;  % regularization factor
    DeltaEps_vec = (A' * A + lambda * eye(Np)) \ (A' * DeltaS_vec);
    % DeltaEps_vec = (A' * A ) \ ( A'*DeltaS_vec);   unstable !!! shit...
    DeltaEps = reshape(DeltaEps_vec, [nx_new, ny_new]);  
 

    % eps_r_predect=eps_r_predect+DeltaEps;  % result to unstability ! er<1 !! shit 
    eps_r_predect=max(real(eps_r_predect+DeltaEps),real(eps_r_bg)) + j*max(imag(eps_r_predect+DeltaEps),0);

    % DeltaEps_itr(itr,:,:)=DeltaEps;  % for debug;
    eps_r_new( Limit_x + 1 : Nx-Limit_x,Limit_x+1 :Nx-Limit_x )= max(eps_r_new(Limit_x+1:Nx-Limit_x ,Limit_x+1 :Nx-Limit_x) + real(DeltaEps),real(eps_r_bg)); %all space
    sigma_new( Limit_x + 1 : Nx-Limit_x,Limit_x+1 :Nx-Limit_x )= max(sigma_new( Limit_x + 1 : Nx-Limit_x,Limit_x+1 :Nx-Limit_x )+ abs(imag(DeltaEps))*eps0,0);

    if( any(eps_r_new(:)<1))
    text='epsr < 1 wtf !!'
    end

    if (mod(itr,5)==0)
        figure(itr+1);
        % subplot(2,1,1);
        % imagesc(real(eps_r_predect)) ;
        % subplot(2,1,2);
        imagesc(real(eps_r_new)) ;
        % caxis([2.846 2.84]]); axis equal tight;
        title(['n=  ' , num2str(itr) ]);
        drawnow;
    end

end
% update er--------------
% Sij_error~= sum of  [ Ei_inc * Ej_inc * delta_er *ds -> find delta er


%% 6. Results and Metrics












%% --- Function Definitions ------------------------------
function [eps_bg,es_inf_cortical, es_inf_normalTrab,eps_true] = define_phantom(Nx,Ny,pml_L,omega)
% Background medium


taw=0.5*1e-12;
eps0 = 8.85e-12

es_inf_bg=2.848; delta_e_bg=1.104; sig_s_bg = 0.005;
es_inf_cortical =8.75; delta_e_cortical=4; sig_s_cortical = 0.01;
es_inf_normalTrab =14; delta_e_normalTrab=7; sig_s_normalTrab = 0.1;

eps_bg = es_inf_bg + delta_e_bg/(1+j*taw*omega) +sig_s_bg/(j*omega*eps0);
eps_crotical= es_inf_cortical + delta_e_cortical/(1+j*taw*omega) +sig_s_cortical/(j*omega*eps0);
eps_normalTrab= es_inf_normalTrab + delta_e_normalTrab/(1+j*taw*omega) +sig_s_normalTrab/(j*omega*eps0);
% True two-layer calcaneus phantom
eps_true = eps_bg * ones(Nx,Ny);
% sig_true = sig_bg * ones(Nx,Ny);

[X,Y] = meshgrid(1:Ny,1:Nx);
center = [Nx/2, Ny/2]; R_out = 20; R_in = 10;
mask = ((X-center(2)).^2+(Y-center(1)).^2) < R_out^2;
eps_true(mask) = eps_crotical;            % cortical
mask2 = ((X-center(2)).^2+(Y-center(1)).^2) < R_in^2;
eps_true(mask2) = eps_normalTrab;           % trabecular
end

function [tx,ty] = circular_array(Nx,Ny,N,pml_L)
theta = linspace(0,2*pi,N+1); theta(end)=[];
% R = min(Nx,Ny)/2-pml_L-50;
R = 35;
cx = Nx/2; cy = Ny/2;
tx = round(cx + R*cos(theta));
ty = round(cy + R*sin(theta));
end

function [E ,Sij ] = fdtd_simulation(Limit_x, nx_new,ny_new,Nx,Ny,dx,dy,Nt,eps_r,sig,k, tx_idx, ty_idx,src,pml_L,sigma_max,sigma_bg,dt ,omega )

% Simplified 2D TM FDTD solver with PEC boundaries
tx=tx_idx(k);
ty=ty_idx(k);
Ntx=size(ty_idx,2);

eps0 = 8.85e-12;
mu0 = 4*pi*1e-7;
c=3e8;

er=zeros(Nx,Ny);
er=eps_r;
% er(cube_i_min :cube_i_max , cube_j_min:cube_j_max )=4;
eps=er*eps0;


%apply PML:
sigma=zeros(Nx,Nt);
sigma=sig;


% PML parameters
order = 3;
sigma_x_1D = zeros(Nx, 1);
sigma_y_1D = zeros(Ny, 1);
kappa_x_1D = ones(Nx, 1);
kappa_y_1D = ones(Ny, 1);

% Grading profile for left and right PML in x-direction
for i = 1:pml_L
    prof = ((pml_L - i + 1)/pml_L)^order;
    sigma_x_1D(i) = sigma_max * prof;
    sigma_x_1D(Nx - i + 1) = sigma_max * prof;
    kappa_x_1D(i) = 1 + (4 - 1) * prof;
    kappa_x_1D(Nx - i + 1) = 1 + (4 - 1) * prof;

    sigma_y_1D(i) = sigma_max * prof;
    sigma_y_1D(Ny - i + 1) = sigma_max * prof;
    kappa_y_1D(i) = 1 + (4 - 1) * prof;
    kappa_y_1D(Ny - i + 1) = 1 + (4 - 1) * prof;
end

[SIGX, SIGY] = meshgrid(sigma_y_1D, sigma_x_1D);
[KAPPAX, KAPPAY] = meshgrid(kappa_y_1D, kappa_x_1D);%%
% sigma=sigma+SIGMA;

% for index=1:pml_L
%     sigma(index,:)=sigma_max*((pml_L-index+2)/pml_L)^1+sigma_bg; %left
%     sigma(Nx-index+1,:)=sigma_max*((pml_L-index+2)/pml_L)^1+sigma_bg; %right
%     sigma(:,index)=sigma_max*((pml_L-index+2)/pml_L)^1+sigma_bg; %up
%     sigma(:,Ny-index+1)=sigma_max*((pml_L-index+2)/pml_L)^1+sigma_bg; %down
% end


%constants for yee algorithem
imagesc(sigma)
R=dt/(2*eps0);
Ra=(c*dt/dx)^2;
Rb=dt/(mu0*dx);

Ca=zeros(Nx,Ny);
Cb=zeros(Nx,Ny);

% For Ezx (x-directed damping)
Ca_x = (1 - R * SIGX ./ er) ./ (1 + R * SIGX ./ er);
Cb_x = Ra ./ ((er + R * SIGX) .* KAPPAX);

% For Ezy (y-directed damping)
Ca_y = (1 - R * SIGY ./ er) ./ (1 + R * SIGY ./ er);
Cb_y = Ra ./ ((er + R * SIGY) .* KAPPAY);
% Ca=(1-R*sigma./eps_r)./(1+R*sigma./eps_r);
% Cb=Ra./(er+R*sigma);

Ez=zeros(Nx,Ny);
Ezx=zeros(Nx,Ny);
Ezy=zeros(Nx,Ny);
Hy=zeros(Nx,Ny);
Hx=zeros(Nx,Ny);

% er draw
figure;
subplot(1,2,1); imagesc(SIGX'); title('\sigma_x') ; axis equal tight;
subplot(1,2,2); imagesc(SIGY'); title('\sigma_y'); axis equal tight;


figure (1);
subplot(2,1,1);
imagesc(real(er));
caxis([1 21]); axis equal tight;
% === Time Loop ===
E=zeros(nx_new,ny_new);
Sij=zeros(1,Ntx);

% v = VideoWriter('fdtd_simulation.avi');  % Name your video file
% v.FrameRate = 20;                        % Set frame rate (frames per second)
% open(v);

for n = 1:Nt
    % --- Update Hx ---
    for i = 1:Nx
        for j = 1:Ny-1
            Hx(i,j) = Hx(i,j) - Rb * (Ez(i,j+1) - Ez(i,j));
        end
    end

    % --- Update Hy ---
    for i = 1:Nx-1
        for j = 1:Ny
            Hy(i,j) = Hy(i,j) + Rb * (Ez(i+1,j) - Ez(i,j));
        end
    end
    % --- Update Ezx ---
    for i = 2:Nx -1
        for j = 1+1:Ny -1
           Ezx(i,j) = Ca_x(i,j)*Ezx(i,j) + (Cb_x(i,j)/Rb)*(Hy(i,j)-Hy(i-1,j));
        end
    end

    % --- Update Ezy ---
    for i = 1+1:Nx -1
        for j = 2:Ny -1

           Ezy(i,j) = Ca_y(i,j)*Ezy(i,j) + (Cb_y(i,j)/Rb)*(Hx(i,j-1)-Hx(i,j));
        end
    end


    Ez=Ezx+Ezy;

    % Ezx(tx,ty) =Ezx(tx,ty) + 10*sin(omega*dt*n);  
    % Ezy(tx,ty) =Ezy(tx,ty) + 10*sin(omega*dt*n);
     Ezx(tx,ty) =Ezx(tx,ty) + src(n) ;  
    % Ezy(tx,ty) =Ezy(tx,ty) + src(n);
    % 

    if( any(isnan(E)))
    text='nan'
    end

    for (antennaIndex=1:size(tx_idx,2))
        Sij(antennaIndex) = Ez(tx_idx(antennaIndex),ty_idx(antennaIndex));   % store Ez field in  the antennas location :D kakkoiiiii
    end
    E=Ez(1+Limit_x:Nx-Limit_x,1+Limit_x:Nx-Limit_x);
    if mod(n,50)==0
        subplot(2,1,2);
        imagesc(real(Ez'));% colorbar;
        title(['n=', num2str(n)]);
        clim([-0.2 0.2]); axis equal tight;

        % draw PML:

        hold on;
        % PML edges (draw rectangles)
        rectangle('Position', [1, 1, Ny, pml_L], 'EdgeColor', 'k', 'LineStyle', '--'); % bottom
        rectangle('Position', [1, Ny-pml_L+1, Ny, pml_L], 'EdgeColor', 'k', 'LineStyle', '--'); % top
        rectangle('Position', [1, 1, pml_L, Nx], 'EdgeColor', 'k', 'LineStyle', '--'); % left
        rectangle('Position', [Ny-pml_L+1, 1, pml_L, Nx], 'EdgeColor', 'k', 'LineStyle', '--'); % right

        %draw circuile 
        xc = Nx/2;     % center x (row index)
        yc = Ny/2;     % center y (column index)
        
        % Draw the circle
        theta = linspace(0, 2*pi, 500);
        x_circle_out = xc + 20 * cos(theta);
        y_circle_out = yc + 20 * sin(theta);
        x_circle_in = xc + 10 * cos(theta);
        y_circle_in = yc + 10 * sin(theta);
        
        % Plot the circle overlay
        plot(y_circle_out, x_circle_out, 'w--', 'LineWidth', 2);  % white dashed circle
         plot(y_circle_in, x_circle_in, 'w--', 'LineWidth', 1);  % white dashed circle
hold off;
        drawnow;
         % frame = getframe(gcf);        % Get current figure as frame
         % writeVideo(v, frame);         % Write frame to video
    end
end
% close (v);

end


function Jacobi = Calc_Jacobi(Np,Ntx,E_inc_k ,ds)

% Np      ---> total number of pixels
% Ntx    ----> Number of transmitters
%A       ----> Jacobi matrix
% E_inc_k  --->  E at k iteration  (E inc : when k=1)


Jacobi = zeros(Ntx^2, Np);  % 256 x Np
row = 1;
for I = 1:Ntx
    Ei = squeeze(E_inc_k(I, :, :));  %   nx_new x ny_new
    for J = 1:Ntx
        Ej = squeeze(E_inc_k(J, :, :));  %  nx_new x ny_new
        product = Ei .* Ej;  % pixel-wise product
        Jacobi(row, :) = reshape(product, [1, Np]);
        row = row + 1;
    end
end

% Multiply by cell area if needed:
Jacobi = Jacobi * ds;

end



