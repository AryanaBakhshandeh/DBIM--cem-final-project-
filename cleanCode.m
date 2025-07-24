
clc;clear all; close all;

% DBIM + FDTD Simulation of Anthropomorphic Two-layer Calcaneus Phantom
%% 1. Domain and Parameters

%physical constants
c=3e8;
eps0 = 8.85e-12;
mu0 = 4*pi*1e-7;


%frequency setting -------------------------------------------
f0 = 1e9;
landa=c/f0;  %waveLength
omega=2*pi*f0;  % radian
T=1/f0;
% grid -------------------------------------------
Nx = 110; Ny = 110;      % grid points

% start_n= 1500;


dx = 1e-3; dy = dx;         %dx = 0.2*landa/(10*sqrt(30)); dy = dx;   % spatial step (1 mm)   dx<landa/(10 sqrt(er_max))  
dt = dx/(sqrt(2*21)*c);         % CFL stability
Nt = round(4*T/dt);               % time steps
Ns=round(Nt) ; % samples of fft

% pml
pml_L=70;
sigma_max=1.5;

%total number of cells
Nxt= Nx+2*pml_L;   Nyt= Nxt;

%antenna and bone radius;
r_bone_in=10 ;% edit kon function ra
r_bone_out=15 ;% edit kon function ra
r_antenna=40;


Nx_box=2*r_bone_out;
 xmin_box=pml_L+0.5*(Nx-Nx_box);
 xmax_box=xmin_box+Nx_box-1



%%source -------------------------------------------

t = (0:Nt-1)*dt;

t0 =0.2e-9;      % Gaussian pulse center and delay
% src = 100*exp(-(t - t0).^2 / (2 * (2e-10)^2)) .*cos(2*pi*f0*t+t0);

T = 0.02e-9;                            % Pulse width [s]

% src = 5* exp(-((t-t0)/T).^2);   % Gaussian pulse

src = 5*sin(omega*t) ;




%% 2. Phantom Definition

[ eps_r_bg,er_inf_cortical, er_inf_normalTrab,eps_r_true] = get_epsPhantom(Nx, Ny,r_bone_in,r_bone_out,omega);

%% 3. Antenna Array Configuration
Ntx = 16;               % transmissions (or 24)
[tx_idx, ty_idx] = place_Antennas(Nxt,Nyt,Ntx,r_antenna);  % antenna's cordinate


% for real bone 
Sij_real=zeros(Ntx,Ntx);

% in forward solber in DBIM
Ek_ij=zeros(Ntx,Nx_box,Nx_box);
Sij=zeros(Ntx,Ntx) ;




%% 4. Precompute Sij in Real case

er_bg_c=real(eps_r_bg);
er_bg=ones(Nx,Ny)*er_bg_c;

sigma_bg_c=abs(imag(eps_r_bg)*omega)*eps0;    % c stands fot constant
sigma_bg=ones(Nx,Ny)*sigma_bg_c;


% Ek_inc = zeros(k,Nx,Ny);
% Sij_inc= zeros (k,k)
sigma_obj=abs(imag(eps_r_true))*omega*eps0;
eps_r_obj=real(eps_r_true);
%incident
% for k = 1:Ntx
%     [E_inc(k,:,:) , Sij_inc(k,:)] =  FDTD(1, Nx,Ny,Nxt,Nyt,dx,dy,Nt,er_bg,sigma_bg,k, tx_idx, ty_idx,src,pml_L,sigma_max,sigma_bg_c,er_bg_c,dt ,omega );
% end
%real bone

 % Sij=zeros(16,16,Ns); %fft
 Sij_real=zeros(Ntx,Ntx);

for k = 1:Ntx
    [trash, Sij_real(k,:)] =  FDTD(Nx,Ny,Nxt,Nyt,dx,dy,Nt,eps_r_obj,sigma_obj,k, tx_idx, ty_idx,src,pml_L,sigma_max,sigma_bg_c,er_bg_c,dt ,omega , r_bone_out,r_bone_in,r_antenna,Ns,Nx_box,xmin_box,xmax_box );
end
clear trash;

% 
% 
%% 5. DBIM + IMATCS Inversion Loop
% %next  : formul 13-2 paynname
% % Sij_predected =Sij_inc+ ?
 Holder_Sij=0;
 ds=dx*dy;
 Holder_Integrand=0;

 er_predect=er_bg_c* ones(Nx_box,Nx_box);
 sig_predect= sigma_bg_c*ones(Nx_box,Nx_box);

 ai=1;
 aj=1;
% % ============================== calculate the Sij scatered : =======================================

% build the Sensivity matrix A-----------------
 eps_r_new=ones(Nx,Ny)*er_bg_c;

 sigma_new =ones(Nx,Ny)*sigma_bg_c;

 eps_r_predect=ones(Nx_box,Nx_box)*eps_r_bg;
 
 Np = Nx_box^2;          % total number of pixels
 A= zeros(Ntx^2, Np);  % 256 x Np

 %start algorithem 

  % er_test=ones(Nx,Ny)*er_bg_c;
  % sigma_test=ones(Nx,Ny)*sigma_bg_c;

for itr=1:40 % iteration

    %first guess
    % if itr>=1
        % forward solver - FDTD
       
        for k = 1:Ntx

            [Ek_ij(k,:,:) , Sij(k,:)] =  FDTD(Nx,Ny,Nxt,Nyt,dx,dy,Nt,eps_r_new,sigma_new,k, tx_idx, ty_idx,src,pml_L,sigma_max,sigma_bg_c,er_bg_c,dt ,omega , r_bone_out,r_bone_in,r_antenna,Ns,Nx_box,xmin_box,xmax_box );
        end
    % end
    % next ->   edit jacobi matrix calculator and ...      cuse fft signal
%     % hase Ns ta  sample !!
    A=Calc_Jacobi(Np,Ntx,Ek_ij,ds);
% 
%     %edit this section !!
    Sij_predected=zeros(Ntx,Ntx);
   
%     %first guess :
    for(I=1:Ntx)
        for(J=1:Ntx)
            Holder_Integrand=0;
            for Ix=1:Nx_box
                for Iy=1:Nx_box
                    contrast = eps_r_predect(Ix,Iy)-(er_bg_c);      % shak daram in bakhsh ra !!!!!!!!!!!!!!!!! $$$$$$$$$$$$$$$$$$$$$ !!! noo noo noo nooo,,,,
                    hold=Ek_ij(I,Ix,Iy).*Ek_ij(J,Ix,Iy);
                    A_size=size(hold);
                    Holder_Integrand=Holder_Integrand+hold*ds *contrast;
                end
            end
            Sij_predected(I,J) = Sij(I,J)-1i*(omega*eps0)/(2*ai*aj).*Holder_Integrand;
        end
    end

    % compute the error
    Sij_error=Sij_real-Sij_predected;
    % sij_holder(itr,:,:)=Sij_predected; % for comparing if we have any change ...

    DeltaS_vec = reshape(Sij_error.', [], 1);  % column vector, size 256x1

    %4. Solve the Inverse Problem (Regularized)
    lambda = 1e-5;  % regularization factor
    DeltaEps_vec = (A' * A + lambda * eye(Np)) \ (A' * DeltaS_vec);
    % DeltaEps_vec = (A' * A ) \ ( A'*DeltaS_vec);   unstable !!! shit...
    DeltaEps = reshape(DeltaEps_vec, [Nx_box, Nx_box]);  


    % eps_r_predect=eps_r_predect+DeltaEps;  % result to unstability ! er<1 !! shit 
    eps_r_predect=max(real(eps_r_predect+DeltaEps),er_bg_c) + j*imag(eps_r_predect+DeltaEps);        % bayad j* imag  ra virayesh konam !!    hes mikonam nadorost ast

    % DeltaEps_itr(itr,:,:)=DeltaEps;  % for debug;
    eps_r_new( xmin_box-pml_L :xmax_box-pml_L ,xmin_box-pml_L :xmax_box-pml_L )= max(eps_r_new(xmin_box-pml_L :xmax_box-pml_L ,xmin_box-pml_L :xmax_box-pml_L) + real(DeltaEps),real(er_bg_c)); %all space
    sigma_new( xmin_box-pml_L :xmax_box-pml_L ,xmin_box-pml_L :xmax_box-pml_L)= max(sigma_new( xmin_box-pml_L :xmax_box-pml_L ,xmin_box-pml_L :xmax_box-pml_L)+ abs(imag(DeltaEps))*eps0*omega,sigma_bg_c);

    if( any(eps_r_new(:)<1))
    text='epsr < 1 wtf !!'
    end

    % if (mod(itr,1)==0 || itr==1)
        figure(itr+1);
        % subplot(2,1,1);
        % imagesc(real(eps_r_predect)) ;
        % subplot(2,1,2);
        imagesc(real(eps_r_new)) ;
        % caxis([2.846 2.84]]); axis equal tight;
        title(['n=  ' , num2str(itr) ]);
        drawnow;
    % end

end
% update er--------------
% Sij_error~= sum of  [ Ei_inc * Ej_inc * delta_er *ds -> find delta er

keyboard;

%% some of functions:

function Jacobi = Calc_Jacobi(Np,Ntx,E_inc_k ,ds)

% Np      ---> total number of pixels in reconstruction region
% Ntx    ----> Number of transmitters
%A       ----> Jacobi matrix
% E_inc_k  --->  E at k iteration  (E inc : when k=1)


Jacobi = zeros(Ntx^2, Np);  % 256 x Np
row = 1;
for I = 1:Ntx
    Ei = squeeze(E_inc_k(I, :, :));  %   Nx_box x Nx_box
    for J = 1:Ntx
        Ej = squeeze(E_inc_k(J, :, :));  %  Nx_box x Nx_box
        product = Ei .* Ej;  % pixel-wise product
        Jacobi(row, :) = reshape(product, [1, Np]);
        row = row + 1;
    end
end

% Multiply by cell area if needed:
Jacobi = Jacobi * ds;

end
% 
