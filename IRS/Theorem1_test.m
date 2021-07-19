%{
Yiliang Liu, School of Cyber Science and Engineering, Xi'an Jiaotong
University, Xi'an 710049, China.

See Theorem 1 in [1].

%} 

clear
close all

nsample = 1e5;

%% system parameters
Ns = 32; Nt = 12; Nr = 4; Ne =1;
SNRt = -30:2:30; sigma = 1; sigma_e = 1; P = 10.^(((SNRt-sigma))./10);
beta = 0.8; alpha = 0.8;
Rs = 3;

%% phase shifter matrix generation
theta = 0+(2*pi).*rand(Ns,1);
Phi = diag(exp(1i.*theta'));

%% main channel estimate and MRT precoding
Gr = sqrt(1/2)*(randn(Nr,Ns)+1i*randn(Nr,Ns));
H = sqrt(1/2)*(randn(Ns,Nt)+1i*randn(Ns,Nt));
Hb = sqrt(1/2)*(randn(Nr,Nt)+1i*randn(Nr,Nt));
[u d v] = svd((alpha*Hb+Gr*Phi*H)'*(alpha*Hb+Gr*Phi*H));
b = u(:,1);

%% wiretap channel sample
He = sqrt(1/2)*(randn(Ne,Nt,nsample)+1i*randn(Ne,Nt,nsample));
Ge = sqrt(1/2)*(randn(Ne,Ns,nsample)+1i*randn(Ne,Ns,nsample));

%% Monte Carlo
for nn = 1:1:nsample
    X(nn) = (beta*He(:,:,nn)*b+Ge(:,:,nn)*Phi*H*b)*(beta*He(:,:,nn)*b+Ge(:,:,nn)*Phi*H*b)';
end

%% calculate secrecy outage probability
for i = 1:length(SNRt)
    Cm = real(log2(1+P(i)*norm((alpha*Hb+Gr*Phi*H)*b)^2));
    x = (2.^(Cm-Rs)-1)/P(i);
    
    % secrecy outage probability expression
    w = beta^2+norm(Phi*H*b)^2;
    p_so(i) = 1-gamcdf(x,Ne,w);
    
    % Monte Carlo for secrecy outage probability
    chi(i) = (2^(Cm-Rs)-1)/P(i);
    p_simu(i) = sum(X>=chi(i))/length(X);
end

p_the1=plot(SNRt,p_so,'r-');
hold on
p_sim1=plot(SNRt,p_simu,'rs');

%% figure set
xlabel('SNR','FontSize',14);
ylabel('Secrecy outage probability','FontSize',14)
leg1= legend({'Theo.','Simu.',},'FontSize',12,'Location','NorthEast');
set(gca,'FontSize',14);



