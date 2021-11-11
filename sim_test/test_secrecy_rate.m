%{
Test AN scheme in 2-2-1 MIMOSE
Liu Yi-Liang
%}
clear; clc; close all
nsample = 1e4;
Nt = 2; Ne = 1; k=1; Nr = 2;
B = 80e6;
snr = 0:50:200;
Eb_N0 = 10.^(((snr))./10);
H = sqrt(k/(2*(k+1)))*(ones(Nr,Nt,nsample)+1i*ones(Nr,Nt,nsample)) + sqrt(1/(2*(k+1)))*(randn(Nr,Nt,nsample)+1i*randn(Nr,Nt,nsample));
he = sqrt(k/(2*(k+1)))*(ones(Ne,1,nsample)+1i*ones(Ne,1,nsample)) + sqrt(1/(2*(k+1)))*(randn(Ne,1,nsample)+1i*randn(Ne,1,nsample));

al = 0.5; be = 1-al;
for j = 1:1:length(Eb_N0)
    for i = 1:1:nsample
        [u d v] = svd(H(:,:,i));
        w = (u(:,1))'; b = v(:,1); z = v(:,2); % w: Bob detection vector; b: Bob beam vector; v: Bob AN precoding
        Cm(i) = log2(1+al*Eb_N0(j)/B*norm(w*H(:,:,i)*b)^2);
        Cw(i) = log2(1+al*Eb_N0(j)/B*norm(he(:,:,i)'*b)^2/(1+be*Eb_N0(j)/B*norm(he(:,:,i)'*z)^2));
        Cs(i) = Cm(i)-Cw(i);
    end
    average_Cs(j) = mean(Cs);
    average_Cm(j) = mean(Cm);
    average_Cw(j) = mean(Cw);
end

throughput = average_Cs(end)*B

figure(1);
subplot(3,1,1);
plot(snr,average_Cs,'r--');
grid on; ylabel('SR (bit/s/Hz)'); xlabel('Eb/N0 (dB)'); title('Ergodic secrecy rate (SR): 80MHz bandwidth, rician k=1, n_t=2, n_r=2, n_e=1');
subplot(3,1,2);
plot(snr,average_Cm,'b--');
grid on; ylabel('BC (bit/s/Hz)'); xlabel('Eb/N0 (dB)'); title('Bob capacity (MC): 80MHz bandwidth, rician k=1, n_t=2, n_r=2');
subplot(3,1,3);
plot(snr,average_Cw,'m--');
grid on; ylabel('EC (bit/s/Hz)'); xlabel('Eb/N0 (dB)'); title('Eve capacity (SR): 80MHz bandwidth, rician k=1, n_t=2, n_e=2');