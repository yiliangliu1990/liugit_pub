function cout =AN_rate(Ne,Nt,s,P)
syms x
f = myeigenpdf(Ne,Nt-s);
n = min(Ne,Nt-s);
cout = 0;
for q = 1:1:n
    cout = cout+double(int(log2(1+P/Nt*x)*f(q),0,Inf));
end
