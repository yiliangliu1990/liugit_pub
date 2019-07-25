%{
Copyright (c) 2019, Yiliang Liu
School of Electronics Information Engineering, Harbin Institute of Technology, China
All rights reserved.
This code is used to provide the Fig. 2 in [1].
[1] XXX

Parameter:
    Ne: the antenna number of Eve.
    Nm: the antenna number of SBS.
%}

clc
clear
close all

syms  k m
nsample  = 1e5;
Pm = 1;
Pv = 1;
Ne = 4;
Nm = 8;
sigma_ve = 1;
sigma_vm = 1;
k_rician = 10; % self-interference channel Rician factor

%% generate channels
h_ve =  sqrt(1/2) * (randn(Ne,1,nsample) + 1i * randn(Ne,1,nsample));
H_mm = sqrt(k_rician/2*(k_rician+1)) * (ones(Nm,Nm,nsample) + 1i * ones(Nm,Nm,nsample))+ sqrt(1/2*(k_rician+1)) * (randn(Nm,Nm,nsample) + 1i * randn(Nm,Nm,nsample));
H_me = sqrt(1/2) * (randn(Ne,Nm,nsample) + 1i * randn(Ne,Nm,nsample));
h_vm =  sqrt(1/2) * (randn(Nm,1,nsample) + 1i * randn(Nm,1,nsample));


%% coding
Cvu = log2(1+Pm*h_vm(:,:,1)'*h_vm(:,:,1)); % setting channel coding based on capacity
Ru = 0:0.1:6; % setting secrecy coding
CCDF_simu = zeros(length(Ru),1);
ccdf_theo = zeros(length(Ru),1);

%% begin simulation
x = Pm/Pv*2.^(Cvu-Ru);
for i = 1:1:nsample
    G_null = null(h_vm(:, :, 1)'*H_mm(:,:,1)); % null-space AN precoding
    z(i) = real(h_ve(:,:,i)'*((sigma_ve/Pm)*eye(Ne)+(H_me(:,:,i)*G_null)*(G_null'*H_me(:,:,i)'))^(-1)*h_ve(:,:,i));
end

for index = 1:length(x)
    th = x(index);
    CCDF_simu(index) = sum(z>=th)/length(z);
end

for i = 1:length(Ru)
    ccdf_theo(i) = exp(-x(i)/(Pm/sigma_ve))*symsum((symsum((factorial(Nm-1)/(factorial(k)*factorial(Nm-1-k)))*x(i)^k,k,0,Ne-m-1)/(1+x(i))^(Nm-1))/(factorial(m))*(x(i)/(Pm/sigma_ve))^m,m,0,Ne-1);
end

%% plot figure
p_mon1 = plot(Ru,CCDF_simu,'rs','LineWidth',1.5,'MarkerIndices',1:4:length(CCDF_simu),'markersize',10);
hold on
p_the1=plot(Ru,ccdf_theo);

