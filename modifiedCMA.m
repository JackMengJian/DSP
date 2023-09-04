function [E_de_x,E_de_y] = modifiedCMA(data_sample_x,data_sample_y,L,mu)

if nargin < 3 %函数输入参数的个数
    L = 9; %L should be odd
end

% Normalize the downsampled input signals  to unit modulus
mx = mean(abs(data_sample_x));
my = mean(abs(data_sample_y));
data_sample_x = data_sample_x/mx;
data_sample_y = data_sample_y/my;
Rxy_in = [data_sample_x;data_sample_y];

% Normalize the input signals to unit modulus
S = size(Rxy_in,2);           
Rxy_temp = rearrange(Rxy_in, L);
gamma = 1;

%Taps initialization
h = zeros(2,2*L);
h(1,round(2*L/4)) = 1;
h(2,round(3*2*L/4)) = 1;
H0 = h;
[ H, e1 ] = depdmcma(Rxy_temp,L,gamma,H0,mu);
hxx = H(1,1:L); hxy = H(2,1:L);

% Once again by CMA
% Update the initial taps
hyy = (hxx').'; hyx = -(hxy').';
H1 = [hxx hyx;hxy hyy];
[ H, e2 ] = depdmcma( Rxy_temp,L,gamma,H1,mu);

Rxy_out = H*Rxy_temp;
E_de_x1 = Rxy_out(1,S-(L-1)/2+1:S);
E_de_y1 = Rxy_out(2,S-(L-1)/2+1:S);
E_de_x = [E_de_x1,Rxy_out(1,1:S-(L-1)/2)];
E_de_y = [E_de_y1,Rxy_out(2,1:S-(L-1)/2)];
save('cma_err.mat','e1','e2');
end