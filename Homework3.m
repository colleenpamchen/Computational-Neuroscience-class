%% Homework3.m 
close all 
clear all


time = [1:101]; 
a=-1;
b=1; 
rnd = a + (b-a).*rand( length(time) ,1); 
gain = 1; % TUNE THIS parameter in the next step  
sig = @(input) 1/( 1 + exp(-1 * gain * input));  

E1 = zeros(length(time) ,1);
E2 = zeros(length(time) ,1);
firing1 = zeros(length(time),1);
firing2 = zeros(length(time),1);

W_extrinsic = 0.1; 
W_intrinsic = 0.1; 
W_inhibitory = 0.1; 
% first iteration, when 
t=1;
N_extrinsic1 = 0.5;
N_extrinsic2 = 0.6;
E1(1)=0;
E2(1)=0;

E1(t+1) = N_extrinsic1 * W_extrinsic + E2(t) * W_intrinsic - E2(t) * W_inhibitory + rnd(t) ; 
firing1(t+1) =  sig(E1(t+1));

E2(t+1) = N_extrinsic2 * W_extrinsic + E1(t) * W_intrinsic - E1(t) * W_inhibitory + rnd(t) ; 
firing2(t+1) = sig(E2(t+1));

for t=2:floor(length(time)/2) % FIRST 50n simulation cycles: 
E1(t+1) = N_extrinsic1 * W_extrinsic + E2(t) * W_intrinsic - E2(t) * W_inhibitory + rnd(t) ; 
firing1(t+1) = firing1(t) + sig(E1(t+1));

E2(t+1) = N_extrinsic2 * W_extrinsic + E1(t) * W_intrinsic - E1(t) * W_inhibitory + rnd(t) ; 
firing2(t+1) = firing2(t) + sig(E2(t+1));

end

% LAST 50 simulation cycles:
W_extrinsic = 2.5;
W_inhibitory = 2.5; 

for t=ceil(length(time)/2):length(time)   
E1(t+1) = N_extrinsic1 * W_extrinsic + E2(t) * W_intrinsic - E2(t) * W_inhibitory + rnd(t) ; 
firing1(t+1) = firing1(t) + sig(E1(t+1));

E2(t+1) = N_extrinsic2 * W_extrinsic + E1(t) * W_intrinsic - E1(t) * W_inhibitory + rnd(t) ;
firing2(t+1) = firing2(t) + sig(E2(t+1));

end

figure
plot(time,firing1(2:102),'k')
hold on 
plot(time,firing2(2:102),'r')
hold off



%%


% keep track of E2 > E1 at 51st 
% and 100th time step 
% proportion of E2>E1 from the 100 runs 
% "activity difference" refers to: EXT2 - EXT1 

%E [:, 2] 


% initialize variables 

% initialize T=1

[]=modulate();


function [E, T50, T100, ratio ] = modulate( E, t, W_intrinsic, W_extrinsic, gain, ActDiff)

% persistent function variables:
N_extrinsic1 = 0.5;
N_extrinsic2 = 0.6;
rnd = a + (b-a).*rand(1);


E(t+1,1) = N_extrinsic1 * W_extrinsic + E2(t) * W_intrinsic - E2(t) * W_inhibitory + rnd  ; % E1

E(t+1,2) = N_extrinsic2 * W_extrinsic + E1(t) * W_intrinsic - E1(t) * W_inhibitory + rnd ; % E2




end













