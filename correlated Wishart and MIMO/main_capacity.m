function cout =main_capacity(Ne,Nt,s,P)
syms x
f = myeigenpdf(Ne,Nt);
cout = 0;
for q = 1:1:s
    cout = cout+double(int(log2(1+P/Nt*x)*f(q),0,Inf));
end
