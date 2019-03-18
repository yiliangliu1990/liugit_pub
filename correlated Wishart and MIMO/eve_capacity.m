function cout =eve_capacity(Ne,Nt,P)
syms x
f = myeigenpdf(Ne,Nt);
n = min(Nt,Ne);
cout = 0;
for q = 1:1:n
    cout = cout+double(int(log2(1+P/Nt*x)*f(q),0,Inf));
end
