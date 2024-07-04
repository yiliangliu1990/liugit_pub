clear
close all

nsample  = 1e2;
SNR = -10:2:10;
Nt = 128;
Nr = 32;
sigma = 1;
Ptot = 10.^(((SNR-sigma))/10);
error = 0.2; % estimated error CN(0,error^2*I)

doa = 30/180*pi; tas = 10/180*pi; dl = 2;
for a = 1:Nt % generate transmit correlation matrix
    for b = 1:Nt
        T(a,b)=exp(-j*2*pi*(b-a)*dl*cos(doa))*exp(-1/2*(2*pi*(b-a)*dl*sin(doa)*tas)^2);
    end
end

for k = 1:length(SNR)
    C_Shannon = 0; C_mmse=0;C_wfl=0;C_zf=0; C_match=0;
    for j = 1:nsample
        H = (sqrt(1/2) * (randn(Nr,Nt) + 1i * randn(Nr,Nt)))*T^(1/2);
        H_estimated = H - error*(sqrt(1/2) * (randn(Nr,Nt) + 1i * randn(Nr,Nt)));
        
        % Shannon
        lambda = eig(H* H');
        pl = WFL(lambda,Ptot(k));
        C_Shannon = C_Shannon+real(sum(log2(1+pl.*lambda)));
        
        % Water-filling
        lambda_e = eig(H_estimated* H_estimated');
        pl_e = WFL(lambda_e,Ptot(k));
        [Wreceive,d,Wtranmit] = svd(H_estimated);
        Pwf = diag(sort(sqrt(pl_e), 'descend'));
        EHwfl = Pwf*Wreceive'*H*Wtranmit;
        for ii = 1:Nr
            Totwfl=EHwfl(ii,:)*(EHwfl(ii,:))';
            Interference_wfl = Totwfl-EHwfl(ii,ii)*EHwfl(ii,ii)';
            C_wfl_se(ii) = real(log2(1+(EHwfl(ii,ii)*EHwfl(ii,ii)')/(1+Interference_wfl)));
        end
        C_wfl =C_wfl+sum(C_wfl_se);
        
        
        % Matching precoding
        Wmatch = H_estimated';%
        for ii = 1:Nr
            Wmatch(:,ii)=Wmatch(:,ii)/norm(Wmatch(:,ii));
        end
        EHmatch=H*Wmatch;
        for ii = 1:Nr
            Totmatch=EHmatch(ii,:)*(EHmatch(ii,:))';
            Interference_match = Totmatch-EHmatch(ii,ii)*EHmatch(ii,ii)';
            C_match_se(ii) = real(log2(1+(Ptot(k)/Nr*EHmatch(ii,ii)*EHmatch(ii,ii)')/(1+Ptot(k)/Nr*Interference_match)));
        end
        C_match =C_match+sum(C_match_se);
        
        % MMSE precoding
        Wmmse = H_estimated'*((H_estimated * H_estimated'+((sigma+error^2*Ptot(k))*Nr/Ptot(k))*eye(Nr)))^(-1);
        for ii = 1:Nr
            Wmmse(:,ii)=Wmmse(:,ii)/norm(Wmmse(:,ii));
        end
        EHmmse=H*Wmmse;
        for ii = 1:Nr
            Tot=EHmmse(ii,:)*(EHmmse(ii,:))';
            Interference = Tot-EHmmse(ii,ii)*EHmmse(ii,ii)';
            C_mmse_se(ii) = real(log2(1+(Ptot(k)/Nr*EHmmse(ii,ii)*EHmmse(ii,ii)')/(1+Ptot(k)/Nr*Interference)));
        end
        C_mmse =C_mmse+sum(C_mmse_se);
        
        % ZF precoding
        Wzf = H_estimated'*(H_estimated* H_estimated')^(-1);
        for ii = 1:Nr
            Wzf(:,ii)=Wzf(:,ii)/norm(Wzf(:,ii));
            %C_zf_se(ii) = log2(1+Ptot(k)/Nr*(abs(H(ii,:)*Wzf(:,ii)))^2);
        end
        EHzf=H*Wzf;
        for ii = 1:Nr
            Tot=EHzf(ii,:)*(EHzf(ii,:))';
            Interference_zf = Tot-EHzf(ii,ii)*EHzf(ii,ii)';
            C_zf_se(ii) = real(log2(1+(Ptot(k)/Nr*EHzf(ii,ii)*EHzf(ii,ii)')/(1+Ptot(k)/Nr*Interference_zf)));
        end
        C_zf = C_zf+sum(C_zf_se);
    end
    average_Shannon(k) = C_Shannon/nsample;
    average_wfl(k) = C_wfl/nsample;
    average_match(k) = C_match/nsample;
    average_mmse(k) = C_mmse/nsample;
    average_zf(k) = C_zf/nsample;
end

plot(SNR,average_Shannon,'k>--');
hold on
plot(SNR,average_wfl,'rd--');
plot(SNR,average_match,'k*--');
plot(SNR,average_mmse,'ms--');
plot(SNR,average_zf,'b^--');
legend({'Shannon capacity','Water-filling','Matching preocding','MMSE precoding','ZF precoding'},'location', 'northwest','FontSize',12);
xlabel('Transmission SNR (dB)');
ylabel('Rate (bit/s/Hz)')
set(gca,'FontSize',16,'FontName','Times New Roman');
set(gca, 'Xlim', [-10 10])
set(gca, 'XTick', -10:2:10)






