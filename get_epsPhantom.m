function [er_bg,es_inf_cortical, es_inf_normalTrab,er_true] = get_epsPhantom(Nxt, Nyt,r_bone_in,r_bone_out,omega)
% Background medium
%parameters 
taw=0.5*1e-12;
eps0 = 8.85e-12

%% eps Parameters _ debay :
% difine debay parameters
es_inf_bg=2.848; delta_e_bg=1.104; sig_s_bg = 0.005;     %background
es_inf_cortical =8.75; delta_e_cortical=4; sig_s_cortical = 0.01;  % bone's cortical  part
es_inf_normalTrab =14; delta_e_normalTrab=7; sig_s_normalTrab = 0.1; % bones' Trab  part  - normal

% calculate the eps_r
er_bg = es_inf_bg + delta_e_bg/(1+1i*taw*omega) +sig_s_bg/(1i*omega*eps0);  
er_crotical= es_inf_cortical + delta_e_cortical/(1+1i*taw*omega) +sig_s_cortical/(1i*omega*eps0);
er_normalTrab= es_inf_normalTrab + delta_e_normalTrab/(1+1i*taw*omega) +sig_s_normalTrab/(1i*omega*eps0);
% True two-layer calcaneus phantom
er_true = er_bg * ones(Nxt,Nyt);



[X,Y] = meshgrid(1:Nyt,1:Nxt);
center = [Nxt/2, Nyt/2];
mask = ((X-center(2)).^2+(Y-center(1)).^2) < r_bone_out^2;
er_true(mask) = er_crotical;            % cortical
mask2 = ((X-center(2)).^2+(Y-center(1)).^2) < r_bone_in^2;
er_true(mask2) = er_normalTrab;           % trabecular
end