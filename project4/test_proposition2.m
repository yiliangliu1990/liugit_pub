%{
Copyright (c) 2019, Yiliang Liu
School of Electronics Information Engineering, Harbin Institute of Technology, China
All rights reserved.
This code is used to provide the Fig. 3 in [1].
[1] XXX

Note:  suitable for high SNR. 

Parameter:
    Ne: the antenna number of Eve.
    Nm: the antenna number of SBS.
%}

clc
clear
close all

syms  k
nsample  = 1e4;
Ne = 4;
Nm = 10;
sigma_mv = 1; % noise in vehicle
sigma_me = 1; % noise in Eve
Cm = zeros(nsample,1); % main capacity
Cw = zeros(nsample,1); % wiretap capacity
Rs = zeros(nsample,1); % instantaneous secrecy rate
average_Rs = zeros(10,1); % average secrecy rate
theoretical_Rs = zeros(10,1); % theoretical ergodic secrecy rate

%% generate channels
H_me = sqrt(1/2) * (randn(Ne,Nm,nsample) + j * randn(Ne,Nm,nsample));
h_mv =  sqrt(1/2) * (randn(1,Nm,nsample) + j * randn(1,Nm,nsample));

%% begin simulation
for n = 1:1:10
    Pm = 10^(((n+20-sigma_mv))/10); % setting transmit SNR = n+20 dBm
    for i = 1:1:nsample
        wd = h_mv(:,:,i)'/norm(h_mv(:,:,i)); % MRT precoding
        Gd = null(h_mv(:,:,i)); % null-space AN precoding
        Cm(i) = log2(1+Pm/Nm*(h_mv(:,:,i)*h_mv(:,:,i)'));
        Cw(i) = log2(det(eye(Ne)+(Pm/Nm*(H_me(:,:,i)*(wd*wd')*H_me(:,:,i)'))/(eye(Ne)+Pm/Nm*(H_me(:,:,i)*(Gd*Gd')*H_me(:,:,i)'))));
        Rs(i) =Cm(i)-Cw(i);
    end
    Rs(Rs<0) = 0; % max(Rs,0)
    average_Rs(n) = real(mean(Rs));
    theoretical_Rs(n) = exp(Nm/Pm)*log2(exp(1))*symsum(expint(k+1,Nm/Pm),k,0,Nm-1)+log2((Nm-Ne)/Nm);
end

SNR_transmit = 20:1:29; % dBm
plot(SNR_transmit, average_Rs,'rs','LineWidth',1.5);
hold on
plot(SNR_transmit, theoretical_Rs,'r-');
