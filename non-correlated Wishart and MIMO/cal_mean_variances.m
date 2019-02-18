%{
Copyright (c) 2016, Yiliang Liu 
Department of Electronics Information Engineering, Harbin Institute of Technology, China
All rights reserved.

This code is used to calculate the contents of Table II in [1].

[1] Secrecy Capacity Analysis of Artificial Noisy MIMO Channels 
- an Approach Based on Ordered Eigenvalues of Wishart Matrices.

Please change Nr and Nt to get the whole contents.

Note that the code is very time-consuming. 
If you have no requirement on accuracy, please use "quad" instead of "int".

Parameter:

    tic: timing the process.
    Nr: the row number of a matrix H.
    Nt: the column number of a matrix H.
    n: the freedom degree of a Whishart matrix HH'(H'H).
    m: the scatter of a Whishart matrix HH'(H'H).
    f(k): the k-th largest eigenvalue pdf of a Whishart matrix HH'(H'H).
    eigen_means(k): the mean of the k-th largest eigenvalue.
    eigen_variances(k): the variance of the k-th largest eigenvalue.

%}

clear all
close all

tic
syms x

Nr = 3;
Nt = 5;
m = max(Nr,Nt);
n = min(Nr,Nt);

eigen_means = zeros(1,n);
eigen_variances = zeros(1,n);

f = myeigenpdf(Nr,Nt); % obtain pdfs of Whishart matrices by myeigenpdf.m

for k = 1:n
    eigen_means(k) = eval(int(f(k)*x,0,Inf));
    eigen_variances(k) = double(simplify(int((x-eigen_means(k))^2*f(k),x,0,inf)));
end

t = toc