%% gozaresh kar :
% changing the Ei and Sij dimentions and delete the time dimention!
%

clc;clear all; close all;
% DBIM + FDTD Simulation of Anthropomorphic Two-layer Calcaneus Phantom

%% 1. Domain and Parameters
Nx = 290; Ny = 290;      % grid points
dx = 1e-3; dy = 1e-3;    % spatial step (1 mm)
Nt = 900;               % time steps
c = 3e8;
dt = dx/(2*c);         % CFL stability


eps0 = 8.85e-12;
mu0 = 4*pi*1e-7;

Nt_sample=Nt/10;

f0 = 1e9; t0 = 1e-9;      % Gaussian pulse center and delay
t = (0:Nt-1)*dt;
source_pulse = exp(-((t-t0)/(t0/4)).^2).*cos(2*pi*f0*t);
w=2*pi*f0;

%pml
pml_L=70;
sigma_max=0.75;

offset=35;
nx_new=Nx-2*(pml_L+offset);
ny_new=Ny-2*(pml_L+offset);

Limit_x=pml_L+offset;

%% 2. Phantom Definition
[ eps_r_bg, sigma_bg, eps_r_true, sigma_true] = define_phantom(Nx, Ny,pml_L);

documentation.define_phantom = @define_phantom;

%% 3. Antenna Array Configuration
Ntx = 16;               % transmissions (or 24)
[tx_idx, ty_idx] = circular_array(Nx, Ny, Ntx,pml_L);

E_inc=zeros(Ntx,nx_new,ny_new);
Sij_true=zeros(Ntx,Ntx);
Sij_inc=zeros(Ntx,Ntx);
documentation.circular_array = @circular_array;

%% 4. Precompute Incident Fields (Background)   Edit this

% E_inc = zeros(Nx, Ny,Nt/10, Ntx);
er_test=ones(Nx,Ny)*eps_r_bg;
sigma_test=ones(Nx,Ny)*sigma_bg;
%   E_inc   %chie neveshtam ?

% Ek_inc = zeros(k,Nx,Ny);
% Sij_inc= zeros (k,k)
for k = 1:Ntx
    [E_inc(k,:,:) , Sij_inc(k,:)] = fdtd_simulation(Limit_x,nx_new,ny_new,Nx, Ny, dx, dy, Nt, er_test, sigma_test,k, tx_idx, ty_idx, source_pulse,pml_L,sigma_max,sigma_bg,dt);
end
documentation.fdtd_simulation = @fdtd_simulation;
% S_true

for k = 1:Ntx
    [ trash , Sij_true(k,:)] = fdtd_simulation(Limit_x,nx_new,ny_new,Nx, Ny, dx, dy, Nt, eps_r_true, sigma_true,k, tx_idx, ty_idx, source_pulse,pml_L,sigma_max,sigma_bg,dt);
end
clear trash;




%% born

%% 5. DBIM + IMATCS Inversion Loop
%next  : formul 13-2 paynname
% Sij_predected =Sij_inc+ ?
Holder_Sij=0;
ds=dx*dy;
Holder_Integrand=0;

eps_r_predect=eps_r_bg * ones(nx_new,ny_new);

ai=1;
aj=1;
% ============================== calculate the Sij scatered : =======================================

% build the Sensivity matrix A-----------------

eps_r_new =er_test;
 Np = nx_new*ny_new;          % total number of pixels
A= zeros(Ntx^2, Np);  % 256 x Np
for itr=1:15 % iteration

    %first guess
    if itr>1  
   % forward solver - FDTD
   er_test=ones(Nx,Ny)*eps_r_bg;
   sigma_test=ones(Nx,Ny)*sigma_bg;
   for k = 1:Ntx
    [E_inc_k(k,:,:) , Sij_inc_k(k,:)] = fdtd_simulation(Limit_x,nx_new,ny_new,Nx, Ny, dx, dy, Nt, eps_r_new, sigma_test,k, tx_idx, ty_idx, source_pulse,pml_L,sigma_max,sigma_bg,dt);
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
    lambda = 1e-5;  % regularization strength
    DeltaEps_vec = (A' * A + lambda * eye(Np)) \ (A' * DeltaS_vec);
    DeltaEps = reshape(DeltaEps_vec, [nx_new, ny_new]);

    eps_r_predect=eps_r_predect+real(DeltaEps);
    % DeltaEps_itr(itr,:,:)=DeltaEps;  % for debug;
    eps_r_new( Limit_x + 1 : Nx-Limit_x,Limit_x+1 :Nx-Limit_x )= eps_r_new(Limit_x+1:Nx-Limit_x ,Limit_x+1 :Nx-Limit_x) + real(DeltaEps); %all space

   
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
function [eps_bg, sig_bg, eps_true, sig_true] = define_phantom(Nx,Ny,pml_L)
% Background medium
eps_bg = 2.848; sig_bg = 0.005;
% True two-layer calcaneus phantom
eps_true = eps_bg * ones(Nx,Ny);
sig_true = sig_bg * ones(Nx,Ny);

[X,Y] = meshgrid(1:Ny,1:Nx);
center = [Nx/2, Ny/2]; R_out = 30; R_in = 15;
mask = ((X-center(2)).^2+(Y-center(1)).^2) < R_out^2;
eps_true(mask) = 12; sig_true(mask) = 0.2;           % cortical
mask2 = ((X-center(2)).^2+(Y-center(1)).^2) < R_in^2;
eps_true(mask2) = 8;  sig_true(mask2) = 0.1;          % trabecular
end

function [tx,ty] = circular_array(Nx,Ny,N,pml_L)
theta = linspace(0,2*pi,N+1); theta(end)=[];
R = min(Nx,Ny)/2-pml_L-4;
cx = Nx/2; cy = Ny/2;
tx = round(cx + R*cos(theta));
ty = round(cy + R*sin(theta));
end

function [E ,Sij ] = fdtd_simulation(Limit_x, nx_new,ny_new,Nx,Ny,dx,dy,Nt,eps_r,sig,k, tx_idx, ty_idx,src,pml_L,sigma_max,sigma_bg,dt )

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

for index=1:pml_L
    sigma(index,:)=sigma_max*((pml_L-index+2)/pml_L)^1+sigma_bg; %left
    sigma(Nx-index+1,:)=sigma_max*((pml_L-index+2)/pml_L)^1+sigma_bg; %right
    sigma(:,index)=sigma_max*((pml_L-index+2)/pml_L)^1+sigma_bg; %up
    sigma(:,Ny-index+1)=sigma_max*((pml_L-index+2)/pml_L)^1+sigma_bg; %down
end
%constants for yee algorithem
R=dt/(2*eps0);
Ra=(c*dt/dx)^2;
Rb=dt/(mu0*dx);

Ca=zeros(Nx,Ny);
Cb=zeros(Nx,Ny);

Ca=(1-R*sigma./eps_r)./(1+R*sigma./eps_r);
Cb=Ra./(er+R*sigma);

Ez=zeros(Nx,Ny);
Ezx=zeros(Nx,Ny);
Ezy=zeros(Nx,Ny);
Hy=zeros(Nx,Ny);
Hx=zeros(Nx,Ny);

% er draw


figure (1);
subplot(2,1,1);
imagesc(real(eps)); axis equal tight;
% === Time Loop ===
E=zeros(nx_new,ny_new);
Sij=zeros(1,Ntx);


for n = 1:Nt
    % --- Update Hx ---
    for i = 2:Nx-1
        for j = 1:Ny-1
            Hx(i,j) = Hx(i,j) - Rb * (Ez(i,j+1) - Ez(i,j));
        end
    end

    % --- Update Hy ---
    for i = 2:Nx-1
        for j = 1:Ny-1
            Hy(i,j) = Hy(i,j) + Rb * (Ez(i+1,j) - Ez(i,j));
        end
    end
    % --- Update Ezx ---
    for i = 2:Nx-1
        for j = 2:Ny-1
            Ezx(i,j) = Ca(i,j)*Ezx(i,j) + (Cb(i,j)/Rb)*(Hy(i,j)-Hy(i-1,j));
        end
    end

    % --- Update Ezy ---
    for i = 2:Nx-1
        for j = 2:Ny-1

            Ezy(i,j) = Ca(i,j)*Ezy(i,j) + (Cb(i,j)/Rb)*(Hx(i,j-1)-Hx(i,j));
        end
    end


    Ez=Ezx+Ezy;

    Ezx(tx,ty) =Ezx(tx,ty) + 10*sin(n*pi/60);
    Ezy(tx,ty) =Ezy(tx,ty) + 10*sin(n*pi/60);


    for (antennaIndex=1:size(tx_idx,2))
        Sij(antennaIndex) = Ez(tx_idx(antennaIndex),ty_idx(antennaIndex));   % store Ez field in  the antennas location :D kakkoiiiii
    end
    E=Ez(1+Limit_x:Nx-Limit_x,1+Limit_x:Nx-Limit_x);
    if mod(n,100)==0
        subplot(2,1,2);
        imagesc(real(Ez'));% colorbar;
        title(['n=', num2str(n)]);
        caxis([-0.5 0.5]); axis equal tight;
        drawnow;
    end
end

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
