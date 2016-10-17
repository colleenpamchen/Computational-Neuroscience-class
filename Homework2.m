%% HOMEWORK 2 
%% Integrate and Fire Neuron 
% V(t+1) = V(t) + (1/Tm)*(El - V(t) + Rm* Ie)

clear all
close all
clc

t = [1:1000]; % ms 
Vreset = -65; % mV
V0= -65; % mV initial condition 
Vth = -50; % THRESHOLD in mV
El = -65; % mV
Tm = 10; % ms
Rm = 10; 
V = NaN(length(t)+1,1); % membrane potential mV 
spikes = NaN(length(t),1);

Ie =[ 1;2;3;6;10];  % array of injected current 
V(1) = -65; % initialized at -65 mV

for ie = 1:length(Ie) 
    
for i = 1:length(t) 
    if V(i) > -50
        V(i+1) = Vreset;
        spikes(i)=1; 
    else
   V(i+1) = V(i)+ ( (1/Tm)*( El - V(i) + Rm * Ie(ie)  )); 
   spikes(i)=0;
    end  
end
spikecount(ie) = sum(nonzeros(spikes)) % frequency: numb of times it spiked 

title('Spike Trains','FontSize',12,'FontWeight','bold')
xlabel('t (ms)','FontSize',12,'FontWeight','bold') 
ylabel('membrane voltage (mV)','FontSize',12,'FontWeight','bold') 

figure;
plot(t(1:100),V(1:100));
legend( strcat('current ', num2str(Ie(ie)), ' mV' ) ,'Location','NorthEast');

end

figure 
plot(Ie, spikecount)
title('FL plot','FontSize',12,'FontWeight','bold')
xlabel('Current','FontSize',12,'FontWeight','bold') 
ylabel('Frequency','FontSize',12,'FontWeight','bold') 


%% Izhikevich Neurons










