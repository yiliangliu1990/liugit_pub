function lowga = Lowgamma(eta,z)
syms t 
lowga = int(t^(eta-1)*exp(-t),0,z);
end

