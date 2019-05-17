%{
chi squared distribution
%}
clear all
close all
syms u
rho = 4;
nsample = 1e4; % Monte Carlo numbers

%% begin Monte Carlo simulations
x =  (randn(rho,1,nsample) + j * randn(rho,1,nsample));
chi = zeros(nsample,1);
for i = 1:1:nsample
    chi(i) = x(:,:,i)'*x(:,:,i);
end

%% theorem expression
f = (u^(rho-1)*exp(-u/2))/(2^(rho)*gamma(rho));

%% figure
h1 = histogram(chi);
h1.Normalization = 'pdf';
[p,ci]=gamfit(chi);
hold on
p_the=ezplot(f,[0,30]);
title('');