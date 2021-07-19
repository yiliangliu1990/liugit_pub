%{
Yiliang Liu, School of Cyber Science and Engineering, Xi'an Jiaotong
University, Xi'an 710049, China.

The CDF of X~X(beta,e,s,u) and its Kolmogorov-Smirnov (KS) test
        x = |beta*a+C*u|^2
        F_X(x) = gamcdf(x,k,w), k = e, w = beta^2+|u|^2

See Lemma 1 in [1].

%} 

clear
close all

nsample = 1e5;

%% distribution parameters
s = 32; e = 5; beta = 0.5; u = randn(s,1);

%% random variables a and C
C = sqrt(1/2)*(randn(e,s,nsample)+1i*randn(e,s,nsample));
a = sqrt(1/2)*(randn(e,1,nsample)+1i*randn(e,1,nsample));

%% generate X via Monte Carlo
for i=1:1:nsample
    X(i) = (beta*a(:,:,i)+C(:,:,i)*u)'*(beta*a(:,:,i)+C(:,:,i)*u);
end

%% Gamma distribution
k = e; w = beta^2+norm(u)^2; x = 0:1:600;
y1 = gamcdf(x,k,w);

%% plot figure
figure;
h1 = histogram(X,'BinWidth',2,'EdgeColor','w');
h1.Normalization = 'cdf';
[p,ci]=gamfit(X);
hold on
plot(x,y1,'r','LineWidth', 2)

%% KS test
y11 = h1.Values;
BinCenters = (h1.BinEdges(1:end-1)+h1.BinEdges(2:end))/2;
y12 = gamcdf(BinCenters,k,w);
max(abs(y12-y11))

%% figure set
xlabel('$x$','FontSize',14,'Interpreter','latex');
ylabel('CDF of $X$','FontSize',14,'Interpreter','latex')
leg1= legend({'CDF of $X\sim X(0,5,32,\mathbf{u})$','Gamma fit of $X\sim X(0,5,32,\mathbf{u})$',},...
    'interpreter','latex','FontSize',12,'Location','SouthEast');
set(gca,'FontSize',14);
set(gca, 'Xlim', [0 600]);
set(gca, 'XTick', 0:100:600);
set(gca, 'Ylim', [0 1]);




