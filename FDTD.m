function [Ek_fft_1ghz_val ,Sk_fft_1ghz_val ] = FDTD(Nx,Ny,Nxt,Nyt,dx,dy,Nt,eps_r,sig,k, tx_idx, ty_idx,src,pml_L,sigma_max,sigma_bg_c,er_bg_c,dt ,omega,rOut,rIn,rAnt,Ns,Nx_box,xmin_box,xmax_box  )

% Simplified 2D TM FDTD solver with PEC boundaries

%  set all Nx    ->   Nxt
%          nx_new --> Nx     ( main grid number of pixels on x axis)

%% constants
eps0 = 8.85e-12;
mu0 = 4*pi*1e-7;
c=3e8;

Ek =zeros(Nx_box,Nx_box,Ns);
Ek_fft=zeros(Nx_box,Nx_box,Ns);
% n_start=Nt+1-Ns; %start sampling;

Skj=zeros(16,Ns); %fft
Skj_time=zeros(16,Ns);

%% define eps and sigma  in all space ( PML + Main grid)
er=ones(Nxt,Nyt)*er_bg_c;
er(pml_L+1:pml_L+Nx, pml_L+1:pml_L+Nx )=eps_r;
eps=er*eps0;

sigma_objct=ones(Nxt,Nyt)*sigma_bg_c;
sigma_objct(pml_L+1:pml_L+Nx, pml_L+1:pml_L+Nx )=sig;

% imagesc(sigma_objct);
% keyboard;

%% antenna's cordination
tx=tx_idx(k);
ty=ty_idx(k);

Ntx=size(ty_idx,2); % number of antennas

%% add PML

m = 3;                             % Polynomial order
% R_0 = 1e-3;                        % Reflection coefficient
% sigma_max = -(m+1)*eps0*c*log(R_0)/(2*pml_L*dx);
%

sigma_x = zeros(Nxt, Nyt);
sigma_y = zeros(Nxt, Nyt);


for i = 1:Nxt
    for j = 1:Nyt

        % Left PML
        if i <= pml_L
            d = (pml_L - i + 1) ;
            sigma_x(i,j) = sigma_max * (d / pml_L)^m +sigma_bg_c;

            % sigma_x(i,j) =3;
        end

        % Right PML
        if i > Nxt - pml_L
            d = (i - (Nxt - pml_L)) ;
            sigma_x(i,j) = sigma_max * (d / pml_L)^m+sigma_bg_c;
            % sigma_x(i,j) =4;
        end

        % Top PML
        if j > Nyt - pml_L
            d = (j - (Nyt - pml_L));
            sigma_y(i,j) = sigma_max * (d / pml_L)^m+sigma_bg_c;
            % sigma_y(i,j)=2;
        end
        % down PML
        if j <= pml_L
            d = (pml_L - j + 1);
            sigma_y(i,j) = sigma_max * (d / pml_L)^m+sigma_bg_c;
            % sigma_y(i,j)=2;
        end
    end
end

% add sigma of bg :
sigma_x=sigma_x+sigma_bg_c;
sigma_y = sigma_y+sigma_bg_c;



% debug sigma -------------------------
% disp('Before debug point');
% figure(1);
% imagesc(sigma_x);
% keyboard; % Enter debug mode here
% figure(2);
% imagesc(sigma_y);
% keyboard;
% disp('After debug point');
% % ------------------------------------------


%% PML Formulation
% Precomputing PML coefficients
be_x = exp(-sigma_x * dt ./eps );
be_y = exp(-sigma_y * dt ./eps);

c_x = zeros(Nxt, Nyt);
c_y = zeros(Nxt, Nyt);

c_mx = zeros(Nxt, Nyt);
c_my = zeros(Nxt, Nyt);




for i = 1:Nxt
    for j = 1:Nyt
        if sigma_x(i,j) > 0
            c_x(i,j) = (1 - be_x(i,j)) / (sigma_x(i,j) * dx);
            c_mx(i,j)=  (1 - be_x(i,j)) / (mu0*sigma_x(i,j) * dx/eps(i,j));
        end
        if sigma_y(i,j) > 0
            c_y(i,j) = (1 - be_y(i,j)) / (sigma_y(i,j) * dy);
            c_my(i,j) = (1 - be_y(i,j)) / (mu0*sigma_y(i,j) * dy/eps(i,j));
        end
    end
end

%for debug --------
% sigma_total = sqrt(sigma_x.^2 + sigma_y.^2);'
%imagesc(sigma_total);
% keyboard;



%% FDTD Time Loop using Yee's scheme

% Initial Fields
Ez = zeros(Nxt, Nyt);
Hx = zeros(Nxt, Nyt);
Hy = zeros(Nxt, Nyt);
Jz = zeros(Nxt, Nyt);

Ezx = zeros(Nxt, Nyt);
Ezy= zeros(Nxt, Nyt);


% Main FDTD loop


isPML=zeros(Nxt,Nyt);
isPML(1:pml_L,:)=1; % left
isPML(Nxt-pml_L+1:Nxt,:)=1;  %right
isPML(:,1:pml_L)=1; isPML(:,Nxt-pml_L+1:Nxt)=1; %down and up


% eta0 = 120*pi;

R=dt/(2*eps0);
Ra=(c*dt/dx)^2;
Rb=dt/(mu0*dx);
Ca = (1-R*sigma_objct./er)./(1+R*sigma_objct./er);
Cb = Ra./(er-R*sigma_objct);

%% ===== Yee..
for n = 1:Nt
    % Update Hx
    for i = 1:Nxt
        for j = 1:Nyt-1
            if(isPML(i,j))
                % dEz_dy = (Ez(i,j+1) - Ez(i,j)) / dy;
                % Hx(i,j) = b_y(i,j) * Hx(i,j) - (dt/mu0) * dEz_dy;
                dEz_y = (Ez(i,j+1) - Ez(i,j)) ;
                Hx(i,j) = be_y(i,j) * Hx(i,j) - (c_my(i,j) * dEz_y);
            else
                Hx(i,j) = Hx(i,j) +Rb*(-1*Ez(i,j+1) + Ez(i,j));
            end

        end
    end

    % Update Hy
    for i = 1:Nxt-1
        for j = 1:Nyt
            if(isPML(i,j))
                dEz_x = (Ez(i+1,j) - Ez(i,j));
                Hy(i,j) = be_x(i,j) * Hy(i,j) + (c_mx(i,j)* dEz_x);
            else
                Hy(i,j) = Hy(i,j) +Rb*(Ez(i+1,j) - Ez(i,j));
            end

        end
    end

    % Adding Source
    % Jz(tx, ty) =  20* src(n);
    % Jz(10, 10) =  20* src(n);

    % Update Ez
    for i = 2:Nxt-1
        for j = 2:Nyt-1
            if(isPML(i,j))
                % dHy_dx = (Hy(i,j) - Hy(i-1,j)) / dx;
                % dHx_dy = (Hx(i,j) - Hx(i,j-1)) / dy;
                % dH = dHy_dx - dHx_dy;
                % Ez(i,j) = be_x(i,j) * be_y(i,j) * Ez(i,j) + (dt / (eps0 * er(i,j))) * (dH - Jz(i,j));
                Ezx(i,j)=be_x(i,j)*Ezx(i,j)+c_x(i,j)*(Hy(i,j)-Hy(i-1,j));
                % if(c_x(i,j) >0)
                % dhy=Hy(i,j)-Hy(i-1,j)
                % end
                Ezy(i,j)=be_y(i,j)*Ezy(i,j)-c_y(i,j)*(Hx(i,j)-Hx(i,j-1));
                % if(c_y(i,j) >0)
                % dhx=(Hx(i,j)-Hx(i,j-1))
                % end
                Ez(i,j) = Ezx(i,j) + Ezy(i,j);
            else

                Ez(i,j) = Ca(i,j).*Ez(i,j) +(1/Rb)*Cb(i,j)*(Hy(i,j)-Hy(i-1,j)  - Hx(i,j)+ Hx(i,j-1) );

            end

        end
    end



    % Source injection
    % Ezx(tx, ty) = Ezx(tx,ty) + src(n);
    % Ezy(tx, ty) = Ezy(tx,ty) + src(n);

    Ez(tx, ty) = Ez(tx,ty) + src(n);

    % E_inc=zeros(Ntx,Nx,Ny);

   

    % store data
    if( n >= Nt+1-Ns)
         %reciveAntenna Data.
        for recIndex=1:16
            recX=tx_idx(recIndex);
            recY=ty_idx(recIndex);
            Skj_time(recIndex,n-(Nt+1-Ns)+1)=Ez(recX,recY);
        end

            % store incident field inside box region
        for Ibox=xmin_box:xmax_box
            for Jbox=xmin_box:xmax_box
                Ek(Ibox-xmin_box+1,Jbox-xmin_box+1,n-(Nt+1-Ns)+1)=Ez(Ibox,Jbox);
            end
        end


    end




    % (Optional) visualize
    if mod(n,1000)==0
        figure(1)
        imagesc(Ez');  title(["Time step: " num2str(n)]);
        caxis([ -0.5 0.5]);  axis equal tight;


        % PML edges (draw rectangles)
        hold on;
        % PML edges (draw rectangles)
        rectangle('Position', [1, 1, Nyt, pml_L], 'EdgeColor', 'k', 'LineStyle', '--'); % bottom
        rectangle('Position', [1, Nyt-pml_L+1, Nyt, pml_L], 'EdgeColor', 'k', 'LineStyle', '--'); % top
        rectangle('Position', [1, 1, pml_L, Nxt], 'EdgeColor', 'k', 'LineStyle', '--'); % left
        rectangle('Position', [Nyt-pml_L+1, 1, pml_L, Nxt], 'EdgeColor', 'k', 'LineStyle', '--'); % right

        %draw circuile
        xc = Nxt/2;     % center x (row index)
        yc = Nyt/2;     % center y (column index)

        % Draw the circle
        theta = linspace(0, 2*pi, 500);
        x_circle_out = xc + rOut * cos(theta);
        y_circle_out = yc + rOut * sin(theta);
        x_circle_in = xc + rIn * cos(theta);
        y_circle_in = yc + rIn* sin(theta);

        % Plot the circle overlay
        plot(y_circle_out, x_circle_out, 'w--', 'LineWidth', 2);  % white dashed circle
        plot(y_circle_in, x_circle_in, 'w--', 'LineWidth', 1);  % white dashed circle
        hold off;
        drawnow;
    end

end

fs=1/dt;
f0=1e9;
df_Ghz=fs/Ns/f0;
n_1ghz=round(1/df_Ghz)+1;

Skj = fft(Skj_time, [], 2)/Ns;  % fft of Skj_time    16xNs
Ek_fft=fft(Ek,[],3)/Ns; % FFT along time (3rd dim), result is complex

Sk_fft_1ghz_val=zeros(1,Ntx);
Ek_fft_1ghz_val=zeros(Nx_box,Nx_box);


% keyboard;

for k=1:16
Sk_fft_1ghz_val(k)=Skj(k,n_1ghz);
end

for I=1:Nx_box
    for J=1:Nx_box
    Ek_fft_1ghz_val(I,J)= Ek_fft(I,J,n_1ghz);
    end
end




end

