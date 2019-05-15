clear all
close all
nsample  = 1e6;
q = randn(nsample,1);
pf = @(x) 1/sqrt(2*pi)*exp(-x^2/2);
cf = @(x) 1/2*(1+erf(x/sqrt(2)));
[c, xxc] = ecdf(q);
[p, xxp]=hist(q);

figure
plot(xxc,c,'r*')
hold on
fplot(cf, [-5 5] ) 
legend({'Monte Carlo', 'Theorem'})

figure
plot(xxp,p/sum(p),'r*')
hold on
fplot(pf, [-5 5] ) 
legend({'Monte Carlo', 'Theorem'})