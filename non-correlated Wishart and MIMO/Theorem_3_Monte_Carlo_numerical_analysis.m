%{
Copyright (c) 2016, Yiliang Liu 
Department of Electronics Information Engineering, Harbin Institute of Technology, China
All rights reserved.

This code is used to provide the secrecy capacity of Theorem 3 in [1].

[1] Secrecy Capacity Analysis of Artificial Noisy MIMO Channels 
- an Approach Based on Ordered Eigenvalues of Wishart Matrices.

Please change Nr and Nt to get corresponding results in Figs.3, 4, 5, 6, 7, 10.
Since the code is very time-consuming, I suggest you record results in
Excel or other texts.

Parameter:

    Nr: the row number of a matrix H.
    Nt: the column number of a matrix H.
    n: the freedom degree of a Whishart matrix HH'(H'H).
    m: the scatter of a Whishart matrix HH'(H'H).
    secrecy_capacity_monte: the Monte Carlo results.
    secrecy_capacity_numerical: the numerical results.
    s: the number of eigen-subchannels for messages. Here, s = n.


%}
clear all;
close all;
syms x l k;
nsample = 1e+04;
Nt = 8;
Nr = 4;
Ne = 4;
m = max(Nr,Nt);
n = min(Nr,Nt);
s = 4; %s = n
sc = zeros(nsample,1);
secrecy_capacity_monte = zeros(nsample,1);
secrecy_capacity_numerical = zeros(nsample,1);
Hm = (randn(Nr,Nt,nsample)+j*randn(Nr,Nt,nsample))/sqrt(2); % main channel
He = (randn(Ne,Nt,nsample)+j*randn(Ne,Nt,nsample))/sqrt(2); % wiretap channel

for p = 1:1:20
    for i = 1:1:nsample;
        P = p*1; % set transmit power
        [u d v] = svd(Hm(:,:,i)'*Hm(:,:,i)); % SVD of main channel
        B = u(:,1:s); % information precoding 
        Z = u(:,(s+1):Nt); % AN signal precoding
        H1=Hm(:,:,i)*B;
        H2=He(:,:,i)*B;
        H3=He(:,:,i)*Z;
        H4=[H2 H3];
        c1(i)=log2(det(eye(Ne)+(P/Nt)*(H4*H4')));
        c2(i)=log2(det(eye(Ne)+(P/Nt)*(H3*H3')));
        c3(i)=log2(det(eye(Nr)+(P/Nt)*(H1*H1')));
        sc(i)=c3(i)+c2(i)-c1(i);
    end
    secrecy_capacity_monte(p) = mean(sc);
    secrecy_capacity_numerical(p) = mycapacity(Nr,Nt,P/Nt)+mycapacity(Ne,Nt-s,P/Nt)-mycapacity(Ne,Nt,P/Nt);
end
x = 1:1:20;
p1 = plot(x,secrecy_capacity_monte(x),'r*','Linewidth',2);
hold on;
p2 = plot(x,secrecy_capacity_numerical(x),'r-','Linewidth',2);
set(gca,'FontSize',14,'FontName','Times New Roman');
set(gca,'xlim',[1 20]);
xlabel('Transmit power (W)','FontSize',16,'FontName','Times New Roman','Interpreter','tex')
ylabel('Average secrecy capacity (bit/Hz/s)','FontSize',16,'FontName','Times New Roman','Interpreter','tex')
legend({'Monte Carlo simulation','Numerical'},'interpreter','latex','FontSize',12,'FontName','Times New Roman');
