%{
Copyright (c) 2017, Yiliang Liu 
Department of Electronics Information Engineering, Harbin Institute of Technology, China
All rights reserved.

The k-th largest eigenvalue distribution (cdf and pdf) of the uncorrelated central Whishart matrix.

This code is used to provide the cdf and pdf in [1].

[1] Secrecy Capacity Analysis of Artificial Noisy MIMO Channels 
- an Approach Based on Ordered Eigenvalues of Wishart Matrices.

Please change Nr and Nt to get the whole contents.

Parameter:

    tic: timing the process.
    Nr: the row number of a matrix H.
    Nt: the column number of a matrix H.
    n: the freedom degree of a Whishart matrix HH'(H'H).
    m: the scatter of a Whishart matrix HH'(H'H).
    cdf(k): the k-th largest eigenvalue cdf of a Whishart matrix HH'(H'H).
    pdf(k): the k-th largest eigenvalue pdf of a Whishart matrix HH'(H'H).
    eigen_means(k): the mean of the k-th largest eigenvalue.
    eigen_variances(k): the variance of the k-th largest eigenvalue.


%}

function yout = myeigenpdf(Nr,Nt)
syms x w
m = max(Nr,Nt);
n = min(Nr,Nt);

s = zeros(n+1,1); % counting the number of different Theta matrices (see [1] Theta matrix definition)
Theta = sym(zeros(n,n,2^n-1)); % 2^n-1 Theta matrices
l = perms(1:1:n);
deter = sym(zeros(1,n));

mu1select1 = zeros(factorial(n),n,n);
mu1select2 = zeros(factorial(n),n,n);
mu2select1 = zeros(factorial(n),n,n);
mu2select2 = zeros(factorial(n),n,n);
selectcom = zeros(factorial(n),n,n);


for k=1:n
    for i = 1:factorial(n)
        mu1select1(i,1:k-1,k) = l(i,1:k-1);
        mu1select2(i,1:k-1,k) = sort(mu1select1(i,1:k-1,k),'ascend');
        mu2select1(i,k:n,k) = l(i,k:n);
        mu2select2(i,k:n,k) = sort(mu2select1(i,k:n,k),'ascend');
    end
end

selectcom = mu1select2+mu2select2;

% matindex is the output mu (see Eqn. (3)) of differet Theta matrices 
for i = 1:n
    matindex{i} = unique(selectcom(:,:,i),'rows'); %unique is used to remove the same items 
end
for k = 1:n
    s(k+1) = s(k)+nchoosek(n,k-1);
end

%{
Output 2^n-1 Theta matrices by subscript i and matindex. 
For example, when k=3, the items of 5th, 6th, and 7th Theta matrices have
column subscript matindex{3}.
%}

for k = 1:n
    for q = 1:nchoosek(n,k-1)
        for i = 1:n
            for j = matindex{k}(q,1:k-1)
                cmd1 = sprintf('igamma(m-n+i+j-1,x)',i,j);
                Theta(i,j,s(k)+q) = eval(cmd1);
            end
            for j = matindex{k}(q,k:n)
                cmd2 = sprintf('Lowgamma(m-n+i+j-1,x)',i,j);
                Theta(i,j,s(k)+q) = eval(cmd2);
            end
        end
    end
end

K=symprod(1/(factorial(m-w)*factorial(n-w)),w,1,n); 

for k = 1:n
    for d = s(1)+1:s(k+1)
        deter(k) = deter(k)+det(Theta(:,:,d));
        cdf(k) = K*deter(k);
        fpdf(k) = diff(cdf(k));   
        yout(k) = fpdf(k);  
    end
end
end
