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

figure
plot(t(1:100),V(1:100));
title('Spike Trains','FontSize',12,'FontWeight','bold')
xlabel('t (ms)','FontSize',12,'FontWeight','bold') 
ylabel('membrane voltage (mV)','FontSize',12,'FontWeight','bold') 
legend( strcat('current ', num2str(Ie(ie)), ' mV' ) ,'Location','NorthEast');

end

figure 
plot(Ie, spikecount)
title('FL plot','FontSize',12,'FontWeight','bold')
xlabel('Current','FontSize',12,'FontWeight','bold') 
ylabel('Frequency','FontSize',12,'FontWeight','bold') 


%% Izhikevich Neurons
% Class 1 Excitable 

clear all
close all
clc


I =[ 1;20;30;60;100];    % I current injection
t = [1:1:500]'; % ms 

v = NaN(length(t)+1,1); % membrane potential mV 
u = NaN(length(t)+1,1);
spikes = NaN(length(t),1);
v(1) =-70;
u(1) =-20;

% class 1 Excitable 
a= 0.02; % a recovery time constasnt
b= -0.1; % b<0 amplifying; b>0 resonating
c=-55; % c spike reset
d= 6; % d outward minus inward currents activated by the spike

figure

for ie = 1:length(I)    
for i = 1:length(t) 
    if v(i) >= 30
       v(i+1)= c;
       u(i+1)= u(i)+d;
       spikes(i)=1; 
    else
       v(i+1) = v(i) + ( 0.04*v(i)^2 + 5*v(i) + 140 - u(i) + I(ie) ); % v membrane potential
       u(i+1) = u(i) + a*( b*v(i) - u(i) ); % u recovery variable
       spikes(i)=0;
    end  
end
spikecount(ie) = sum(nonzeros(spikes)) % frequency: numb of times it spiked 

title('Spike Trains for Class 1 Excitable','FontSize',12,'FontWeight','bold')
xlabel('t (ms)','FontSize',12,'FontWeight','bold') 
ylabel('membrane voltage (mV)','FontSize',12,'FontWeight','bold') 

subplot(length(I),1,ie)
plot(t(1:100),v(1:100));
legend( strcat('current ', num2str(I(ie)), ' mV' ) ,'Location','NorthEast');
end

figure 
plot(I, spikecount)
title('FL plot for Class 1 Excitable','FontSize',12,'FontWeight','bold')
xlabel('Current','FontSize',12,'FontWeight','bold') 
ylabel('Frequency','FontSize',12,'FontWeight','bold') 


%% Class 2 Excitable
clear all
close all
clc
% class 2 Excitable 
a= 0.2;% a recovery time constasnt
b= 0.26;% b<0 amplifying; b>0 resonating
c= -65; % c spike reset
d= 0;% d outward minus inward currents activated by the spike

I =[ 1;20;30;60;100];    % I current injection
t = [1:1:500]'; % ms 

v = NaN(length(t)+1,1); % membrane potential mV 
u = NaN(length(t)+1,1);
spikes = NaN(length(t),1);
v(1) =-70;
u(1) =-20;

figure

for ie = 1:length(I)    
for i = 1:length(t) 
    if v(i) >= 30
       v(i+1)= c;
       u(i+1)= u(i)+d;
       spikes(i)=1; 
    else
       v(i+1) = v(i) + ( 0.04*v(i)^2 + 5*v(i) + 140 - u(i) + I(ie) ); % v membrane potential
       u(i+1) = u(i) + a*( b*v(i) - u(i) ); % u recovery variable
       spikes(i)=0;
    end  
end
spikecount(ie) = sum(nonzeros(spikes)) % frequency: numb of times it spiked 

title('Spike Trains for Class 2 Excitable','FontSize',12,'FontWeight','bold')
xlabel('t (ms)','FontSize',12,'FontWeight','bold') 
ylabel('membrane voltage (mV)','FontSize',12,'FontWeight','bold') 

subplot(length(I),1,ie)
plot(t(1:100),v(1:100));
legend( strcat('current ', num2str(I(ie)), ' mV' ) ,'Location','NorthEast');
end

figure 
plot(I, spikecount)
title('FL plot for Class 2 Excitable','FontSize',12,'FontWeight','bold')
xlabel('Current','FontSize',12,'FontWeight','bold') 
ylabel('Frequency','FontSize',12,'FontWeight','bold') 

%% Regular spiking neuron 
clear all
close all
clc

a= 0.02;
b= 0.2;
c= -65;
d= 8;

I =[ 1;20;30;60;100];    % I current injection
t = [1:1:500]'; % ms 

v = NaN(length(t)+1,1); % membrane potential mV 
u = NaN(length(t)+1,1);
spikes = NaN(length(t),1);
v(1) =-70;
u(1) =-20;

figure

for ie = 1:length(I)    
for i = 1:length(t) 
    if v(i) >= 30
       v(i+1)= c;
       u(i+1)= u(i)+d;
       spikes(i)=1; 
    else
       v(i+1) = v(i) + ( 0.04*v(i)^2 + 5*v(i) + 140 - u(i) + I(ie) ); % v membrane potential
       u(i+1) = u(i) + a*( b*v(i) - u(i) ); % u recovery variable
       spikes(i)=0;
    end  
end
spikecount(ie) = sum(nonzeros(spikes)) % frequency: numb of times it spiked 

title('Spike Trains for Regular spiking neuron','FontSize',12,'FontWeight','bold')
xlabel('t (ms)','FontSize',12,'FontWeight','bold') 
ylabel('membrane voltage (mV)','FontSize',12,'FontWeight','bold') 

subplot(length(I),1,ie)
plot(t(1:100),v(1:100));
legend( strcat('current ', num2str(I(ie)), ' mV' ) ,'Location','NorthEast');
end

figure 
plot(I, spikecount)
title('FL plot for Regular spiking neuron','FontSize',12,'FontWeight','bold')
xlabel('Current','FontSize',12,'FontWeight','bold') 
ylabel('Frequency','FontSize',12,'FontWeight','bold') 









