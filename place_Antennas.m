
function [tx,ty] = place_Antennas(Nxt,Nyt,Ntx,r_antenna)
theta = linspace(0,2*pi,Ntx+1); theta(end)=[];
cx = Nxt/2; cy = Nyt/2;

tx = round(cx + r_antenna*cos(theta));
ty = round(cy + r_antenna*sin(theta));
end
