function cout = myeigencapacity(Nr,Nt,k,p)
syms x 
f = myeigenpdf(Nr,Nt); % obtain pdfs of Whishart matrices by myeigenpdf.m
cout = double(int(log2(1+p*x)*f(k),0,Inf));
end