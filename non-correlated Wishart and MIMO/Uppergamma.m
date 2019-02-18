function upperga = Uppergamma(eta,z)
syms t 
upperga = int(t^(eta-1)*exp(-t),z,Inf);
end

