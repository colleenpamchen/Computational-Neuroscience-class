%% FINAL EXAM: Colleen Chen 
%  
% [ ]  FIRST try to make the STDP + Homeostasis learning rule work with
%       non random stimulus.  
% ( ) alpha- homeostatic scaling factor??  
% ( ) R average firing rate of the postsynpatic neuron j
% ( ) gamma- tuning factor =50
%

clear all
close all
clc

% ADD NEW homeostatic variables: 
alpha = 10; % homeostatic scaling factor? (not given) 
R = 0; % average firing rate of postsynaptic neuron j = OUTPUT neuron 
T = 5000 ; % ms. %T=5 SECONDS. time scale over which the firing rate of the postsunaptic  neuron was averaged 
Rtarget= 35; % 35 Hz target firing rate 
gamma = 50; % homeostasis tuning factor 

% ( alpha * S1(i)* (1-(R1/ Rtarget) ) + ( LTP1(i)+LTD1(i) )  )

Ne=100;       Ni=0;     Nout=1;
N = Ne + Ni;

% parameters of the excitatory RS Izzy neuron 
a= 0.02;
b= 0.2;
c= -65;
d= 8;

% S contains the weights ordered from, to
% initialize weights according to network1 configuration 
% IMPLEMENT A all-to-one connections 
S1 = zeros(N,1); %% ?????????????? for the weight matrix, is it a
% [poisson neurons 100 x 1 izzy output neuron] ????
% becaues if [N N] then it implies all-to-all connection, which is 
% how izzy had it set up, and later, in the 
aa=.01; bb=.03;
for i = 1:N
    for j = 1:Nout
        S1(i,j) = (aa + (bb-aa).*rand(1,1) ) ; %.* 1000;  
    end
end

% STDP lists for LTD and LTP weight changes
LTP1(1:Ne, 1:Nout) = 0;    % ordered from, to
LTD1(1:Ne, 1:Nout) = 0;    % ordered from, to

% STDP parameters
A_plus = 0.0002;
A_minus = 0.000066;
t_plus = 20;
t_minus = 60;

wmin = 0;
wmax = 0.03;

Time=1000;  % conveninetly, both simulation seconds and timestep ms 

% there is only 1 v 
v1 = -65 ;
u1 = b*v1;  % Initial values of u

fr1 = zeros(Time,1); % FIRING RATE keeps track of how many cycles per second output neuron fired 
spikes1 = zeros(N,Time); % input POISSON spikes into the system 
I1=zeros(N,1); 

% keep track of the time steps when the poisson spiked. 
lambda=[0.2:0.2:20]; % mean Firing Rate of poisson neurons 
% this would be the value of pixels, all 28*28 values. W= (64) izzy neuron
% would have 784 weights [784, 64]
% H= firing rate of 64 neurons  [64, 1] 
% [784,64] * [64,1] = [784, 1] image 
for i=1:N
   for t=1:1000
      xrand = rand(1);
      % this is generating 'tau' the interspike intervals 
      spikes1(i,t+1) =  spikes1(i,t) - ( log(xrand)/lambda(i) )*1000 ;
   end
end
spikes1 = ceil(spikes1);

for sec = 1:Time  % 1000 simulation seconds 
    disp(sec);
    vfired1=[]; 
    % this gets reset at every second loop 
    for t=1:1000          % simulation of 1000 ms
        
        for i=1:N
                if v1 >=30 
                    fr1(sec)=fr1(sec)+1;
                    v1 = c;
                    u1 = u1 + d;
                    vfired1 = [vfired1; t]; % v neuron fired at t 
                
                S1(i) = max(wmin, S1(i) + LTD1(i)); % LTD to weights that connect to all exc from i
                LTD1(i) = A_minus;  % set max LTD to i from all exc
                end
                
                Ifired1 = find(spikes1(i,1:Time+1)==t); % the t indices of ith poisson that fired 
                if ~isempty(Ifired1) 
                        I1(i) = I1(i) + S1(i);  %*length(fired) ; % If at this timestep, neuron i spiked, then add input and synaptic weight
                       
                        S1(i) = min(wmax, S1(i)+ LTP1(i)); % LTP to weights that connect to i from all exc 
                        LTP1(i) = A_plus;   % set max LTP to all exc from i          
                end
        v1 = v1 + 0.5*(0.04* v1^2 +5*v1 +140 -u1 + I1(i) );
        v1 = v1 + 0.5*(0.04* v1^2 +5*v1 +140 -u1 + I1(i) );
        u1 = u1 + a.*(b.*v1-u1);
    
        end % i neuron loop 
         % exponentially decay LTD and LTP based on time constants
        LTP1 = LTP1 - LTP1/t_plus;
        LTD1 = LTD1 - LTD1/t_minus;
        
    end % t ms loop
    
end % loop: Time second 

% figure
subplot(2,2,1)
plot((1:Time), fr1./T) % time, firing rate
title('HOMEOSTASIS')
% figure
subplot(2,2,2)
plot((1:100),S1,'+r')    % synapse ID number, synaptic strength 
title('HOMEOSTASIS')
ylim([-0.01 0.06])    


% Standard Izzy STDP 
% FIRST NETWORK CONFIG WITHOUT HOMEOSTASIS: 
% IMPLEMENT A all-to-one connections: 

S = zeros(N,1); % [poisson neurons 100 x 1 izzy output neuron] 
aa=.01; bb=.03;
for i = 1:N
        S(i) = (aa + (bb-aa).*rand(1,1) ) ; % Synaptic weights 100 x 1 outputNeuron 
end

% STDP lists for LTD and LTP weight changes
LTP(1:Ne, 1:Nout) = 0;    % ordered from, to
LTD(1:Ne, 1:Nout) = 0;    % ordered from, to

% there is only 1 v 
v = -65; 
u = b*v;  % Initial values of u

fr = zeros(Time,1); % FIRING RATE keeps track of how many cycles per second output neuron fired 
spikes = zeros(N,Time); % input POISSON spikes into the system 
fired=[];
I=zeros(N,1); 
% nspikes = zeros(N,1);

lambda=[0.2:0.2:20]; % mean Firing Rate of poisson neurons 
for i=1:N
   for t=1:1000 % I guess this is looping through the Time trials from 1:1000 SECONDS
      xrand = rand(1);
      % this is generating 'tau' the interspike intervals 
      spikes(i,t+1) =  spikes(i,t) - ( log(xrand)/lambda(i) )*1000 ; % generate 1000 spikes per neuron in SECONDs  
   end
end
spikes = ceil(spikes); % these are the spike times in ms.  

for sec = 1:Time  % 1000 simulation seconds 
    disp(sec);
    vfired=[]; 
    
    for t=1:1000 %1000          % simulation of 1000 ms
       
        for i = 1:N  % this loops through the poisson input NEURONS 
          
            if v >=30
                fr(sec)=fr(sec)+1;
                v = c;
                u = u + d;   
                vfired = [vfired; t]; % v neuron fired at t 
                
                S(i) = max(wmin, S(i) + LTD(i)); % LTD to weights that connect to all exc from i
                LTD(i) = A_minus;  % set max LTD to i from all exc
            end 
            
                Ifired = find(spikes(i,1:Time+1)==t); % the t indices of ith poisson that fired 
                if ~isempty(Ifired) 
                        I(i) = I(i) + S(i);  %*length(fired) ; % If at this timestep, neuron i spiked, then add input and synaptic weight
                       
                        S(i) = min(wmax, S(i)+ LTP(i)); % LTP to weights that connect to i from all exc 
                        LTP(i) = A_plus;   % set max LTP to all exc from i          
                end
             v = v + 0.5* (0.04 * v^2 + 5*v + 140 -u + I(i) );
             v = v + 0.5*(0.04 * v^2 + 5*v + 140 -u + I(i) );       
             u = u + a*(b*v -u ); 
        end % poissonNeuron i loop
        LTP = LTP - LTP/t_plus;
        LTD = LTD - LTD/t_minus;
    end % end of t ms loop    
    
end % end of Time s loop 

% figure
subplot(2,2,3)
plot([1:Time], fr./5 ) % time, firing rate
% figure
subplot(2,2,4)
plot([1:100],S,'+r')    % synapse ID number, synaptic strength 
ylim([-0.01 0.06])    

%%
% Neuron network setup:
%   100 poisson spiking neurons
%   Excitatory Izhekevich regular spiking (RS)
%   Inhibitory fast spiking (FS)
%   delay: 1ms 
%   initial weights: randn()/uniform distribution [.01 .03]
%   STDP paramters: 
%       apre, apost, taupre, taupost 
% 
% 
% synaptic weights
% firing rates
% target homeostatsis firing rate: 35hz 
% 
% Homeostasis & STDP Learning rule:   
% 
% weight update
% 
% 
% 
% 
% poisson_spikes = zeros(N,1000); 
% poisson_spikes(1,1)=0;
% for n=1:100
%    for t=1:1000-1
%        xrand = rand(1);
%       poisson_spikes(n,t+1) =  poisson_spikes(n,t) - log(xrand)/lambda(i) ;
% %       fired = find(poisson_spikes>xrand); 
%    end
% end
% 
%  
% % 
% % 
% % % 
% for i = 1:N
% poisson_spikes(i,:) = poissrnd( lambda(i),1,1000 );
% end
% 
% % 
% 
