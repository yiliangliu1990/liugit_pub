function pl = WFL(lambda, Ptot)   
[lambda, idx] = sort(lambda, 'descend');
lambda = lambda(lambda > 0);      
pl = -1;
try
    while (min(pl) < 0)
        mu = (Ptot + sum(1 ./ lambda)) / length(lambda);
        pl = mu - 1 ./ lambda;
        lambda = lambda(1:end-1);
    end
catch
    disp('There exists no water filling level for the input eigenvalues. Check your data and try again')
end
pl = [pl; zeros(length(idx) - length(pl), 1)]; % assigning zero power for weak eigen-modes 
pl(idx) = pl; % rearranging the power levels 