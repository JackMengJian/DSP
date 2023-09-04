function [X_out,Y_out] = NLC_PDM(X_in,Y_in,gama,L)
 leff = (1-exp(-2*0.2e-3/4.342945*L))/(2*0.2e-3/4.342945);
% leff = L;
r = gama*1e-3; 
phase = 1*(8/9*r*(abs(X_in).^2+abs(Y_in).^2)*leff);
X_out = exp(1i*phase).*X_in;
Y_out = exp(1i*phase).*Y_in;
end

