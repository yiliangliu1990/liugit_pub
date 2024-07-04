clear
close all

nsample  = 1e1;
SNR = -10:2:10;
Nr = 128; %基站天线
K = 32; % 用户
sigma = 1;
Ptot = 10.^(((SNR-sigma))/10);


for k = 1:length(SNR)
    C_mmse=0;C_zf=0; C_match=0;
    for j = 1:nsample
        H = (sqrt(1/2) * (randn(Nr,K) + 1i * randn(Nr,K)));
       
        % Matching precoding
        Wmatch = H';% K X Nr
        for ii = 1:K
            Wmatch(ii,:)=Wmatch(ii,:)/norm(Wmatch(ii,:));
        end
        EHmatch=Wmatch*H;
        for ii = 1:K
            Totmatch=(EHmatch(:,ii))'*EHmatch(:,ii);
            Interference_match = Totmatch-EHmatch(ii,ii)*EHmatch(ii,ii)';
            C_match_se(ii) = real(log2(1+(Ptot(k)/K*EHmatch(ii,ii)*EHmatch(ii,ii)')/(1+Ptot(k)/K*Interference_match)));
        end
        C_match =C_match+sum(C_match_se);
        
        % MMSE precoding
        Wmmse = (H' * H+(sigma+(K/Ptot(k))*eye(K)))^(-1)*H'; 
        for ii = 1:K
            Wmmse(ii,:)=Wmmse(ii,:)/norm(Wmmse(ii,:));
        end
        EHmmse=Wmmse*H;
        for ii = 1:K
            Tot=(EHmmse(:,ii))'*EHmmse(:,ii);
            Interference = Tot-EHmmse(ii,ii)*EHmmse(ii,ii)';
            C_mmse_se(ii) = real(log2(1+(Ptot(k)/K*EHmmse(ii,ii)*EHmmse(ii,ii)')/(1+Ptot(k)/K*Interference)));
        end
        C_mmse =C_mmse+sum(C_mmse_se);
        
        % ZF precoding
        Wzf = (H'*H)^(-1)*H';
        for ii = 1:K
            Wzf(ii,:)=Wzf(ii,:)/norm(Wzf(ii,:));
            %C_zf_se(ii) = log2(1+Ptot(k)/Nr*(abs(H(ii,:)*Wzf(:,ii)))^2);
        end
        EHzf=Wzf*H;
        for ii = 1:K
            Tot=(EHzf(:,ii))'*EHzf(:,ii);
            Interference_zf = Tot-EHzf(ii,ii)*EHzf(ii,ii)';
            C_zf_se(ii) = real(log2(1+(Ptot(k)/K*EHzf(ii,ii)*EHzf(ii,ii)')/(1+Ptot(k)/K*Interference_zf)));
        end
        C_zf = C_zf+sum(C_zf_se);
    end
    average_match(k) = C_match/nsample;
    average_mmse(k) = C_mmse/nsample;
    average_zf(k) = C_zf/nsample;
end




plot(SNR,average_match,'k*--');
hold on
plot(SNR,average_mmse,'ms--');
plot(SNR,average_zf,'b^--');
legend({'Matching preocding','MMSE precoding','ZF precoding'},'location', 'northwest','FontSize',12);
xlabel('Transmission SNR (dB)');
ylabel('Rate (bit/s/Hz)')
set(gca,'FontSize',16,'FontName','Times New Roman');
set(gca, 'Xlim', [-10 10])
set(gca, 'XTick', -10:2:10)






