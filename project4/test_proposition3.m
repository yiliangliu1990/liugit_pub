%{
Copyright (c) 2019, Yiliang Liu
School of Electronics Information Engineering, Harbin Institute of Technology, China
All rights reserved.
This code is used to provide the offloading simulations in [1].
[1] XXX
Parameter:
%}

clear
close all

%% vehicle computing parameters
av = 5; % computing capacity of vehicle 5 bit/s/Hz
vehicle_core = 2e9; % CPU of vehicle 1*2 Ghz
vehicle_computing_capacity = av*vehicle_core;

%% MEC computing parameters
am = 1:0.5:4; % computing capacity of MEC bit/s/Hz
MEC_core = 64e9; % CPU of MEC 8*64Ghz
MEC_computing_capacity = am*MEC_core;

%% communication throughput
Ru = 5; % uplink secrecy rate bit/s/Hz
Rd = 8; % downlink secrecy rate bit/s/Hz
Rd_sdp = 12; % downlink secrecy rate with SDP bit/s/Hz
bandwidth_uplink = 5e6; % 50 MHz
bandwidth_downlink = 2e6; % 10 MHz
modulation = 1024; % 1024 OFDM modulation
uplink_throughput = Ru*bandwidth_uplink*modulation;
downlink_throughput = Rd*bandwidth_downlink*modulation;
downlink_throughput_sdp = Rd_sdp*bandwidth_downlink*modulation;

%% data
alpha = 0.8; % ratio of input and output of MEC
m  = 128; % 1 Mb
package_size = m*8*1024*1024; 

T1 = zeros(length(am),1);
T1_sdp = zeros(length(am),1);
Tv = zeros(length(am),1);
Tm = zeros(length(am),1);

%% begin simulation
for i = 1:1:length(am)
    a1 = 1/MEC_computing_capacity(i)+1/uplink_throughput+alpha/downlink_throughput;
    a1_sdp = 1/MEC_computing_capacity(i)+1/uplink_throughput+alpha/downlink_throughput_sdp;
    eta = 1-1/(vehicle_computing_capacity*a1+1);
    eta_sdp = 1-1/(vehicle_computing_capacity*a1_sdp+1);
    T1(i) = package_size*eta/vehicle_computing_capacity;
    T1_sdp(i) = package_size*eta_sdp/vehicle_computing_capacity;
    Tv(i) = package_size/vehicle_computing_capacity;
    Tm(i)= package_size/MEC_computing_capacity(i)+package_size/uplink_throughput+package_size*alpha/downlink_throughput;
end

plot(am,T1*1000,'rd-','Linewidth',1.5,'markersize',10)
hold on
plot(am,Tv*1000,'b^-','Linewidth',1.5,'markersize',10)
plot(am,Tm*1000,'ms-','Linewidth',1.5,'markersize',10)
plot(am,T1_sdp*1000,'k--','Linewidth',1.5,'markersize',10)

