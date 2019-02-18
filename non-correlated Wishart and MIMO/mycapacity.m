function cout = mycapacity(Nr,Nt,p)
syms x k v
m = max(Nt,Nr);
n = min(Nt,Nr);
cout = double(int(log2(1+p*x)*symsum((factorial(k)/factorial(k+m-n))*(symsum((-1)^v*nchoosek(k+m-n,k-v)*(x^v/factorial(v)),v,0,k))^2,k,0,n-1)*exp(-x)*x^(m-n),0,Inf));
end