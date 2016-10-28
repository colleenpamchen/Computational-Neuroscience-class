%% Homework3.m 
close all 
clear all


time = [1:101]; 
a=-1;
b=1; 
rnd = a + (b-a).*rand( length(time) ,1); 
gain = 1; 
sig = @(input) 1/( 1 + exp(-1 * gain * input));  

E1 = zeros(length(time) ,1);
E2 = zeros(length(time) ,1);

W_extrinsic = 0.1; 
W_intrinsic = 0.1; 
W_inhibitory = 0.1; 
% first iteration, when 
t=1;
N_extrinsic1 = 0.5;
N_extrinsic2 = 0.6;
E1(1)=1;
E2(1)=1;

E1(t+1) = N_extrinsic1 * W_extrinsic + E2(t) * W_intrinsic - E2(t) * W_inhibitory + rnd(t) ; 
E2(t+1) = N_extrinsic2 * W_extrinsic + E1(t) * W_intrinsic - E1(t) * W_inhibitory + rnd(t) ; 

for t=2:floor(length(time)/2) % FIRST 50n simulation cycles: 
E1(t+1) = N_extrinsic1 * W_extrinsic + sig(E2(t)) * W_intrinsic - sig(E2(t)) * W_inhibitory + rnd(t) ; 
E2(t+1) = N_extrinsic2 * W_extrinsic + sig(E1(t)) * W_intrinsic - sig(E1(t)) * W_inhibitory + rnd(t) ; 
end

% LAST 50 simulation cycles:
W_extrinsic = 2.5;
W_inhibitory = 2.5; 

for t=ceil(length(time)/2):length(time)  
E1(t+1) = N_extrinsic1 * W_extrinsic + sig(E2(t)) * W_intrinsic - sig(E2(t)) * W_inhibitory + rnd(t) ; 
E2(t+1) = N_extrinsic2 * W_extrinsic + sig(E1(t)) * W_intrinsic - sig(E1(t)) * W_inhibitory + rnd(t) ; 
end

figure
plot(time,E1(2:102),'k')
hold on 
plot(time,E2(2:102),'r')
hold off



