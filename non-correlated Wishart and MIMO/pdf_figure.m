%{
Copyright (c) 2016, Yiliang Liu 
Department of Electronics Information Engineering, Harbin Institute of Technology, China
All rights reserved.

This code is used to provide Fig. 1 in [1].

[1] Secrecy Capacity Analysis of Artificial Noisy MIMO Channels 
- an Approach Based on Ordered Eigenvalues of Wishart Matrices.

Please change Nr and Nt to get other pdf figures of uncorrelated central Whishart matrices.

Parameter:

    Nr: the row number of a matrix H.
    Nt: the column number of a matrix H.
    n: the freedom degree of a Whishart matrix HH'(H'H).
    m: the scatter of a Whishart matrix HH'(H'H).
    f(k): the k-th largest eigenvalue pdf of a Whishart matrix HH'(H'H).


%}
clear all
close all


syms x

Nr = 4;
Nt = 7;
m = max(Nr,Nt);
n = min(Nr,Nt);

f = myeigenpdf(Nr,Nt); % obtain pdfs of Whishart matrices by myeigenpdf.m

for i = 1:n
    p=ezplot(f(i),[0,25,0,0.7]);
    set(p,'linewidth',2);
    hold all;
end
title('');