%% HOMEWORK 2 
%% Integrate and Fire Neuron 
% V(t+1) = V(t) + (1/Tm)*(El - V(t) + Rm* Ie)

clear all
close all
clc

t = [1:1000]; % ms 
Vreset = -65; % mV
V0= -65; % mV initial condition 
Vth = -50; % mV
El = -65; % mV
Tm = 10; % ms
Rm = 10; 
V = NaN(length(t)+1,1); 
spikecount = NaN(length(t),1);

Ie =[ 1.6 ]; %;20;30;40;50];  % array of injected current 

% for V0 initialized at -65 mV
V(1) = -65;  

for i = 1:length(t) 
   
    if V(i) > -50
        V(i+1) = Vreset;
        spikecount(i)=1; 
    else
   V(i+1) = V(i)+ ((1/Tm)*( El - V(i) + Rm*Ie)); 
   spikecount(i)=0;
    end
    
end

plot(t(1:100),V(1:100))

sum(nonzeros(spikecount)) % frequency: numb of times it spiked 

%% Izhikevich Neuron 
