%% FINAL EXAM: Colleen Chen 
%  
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
% 
% 
% 
%% Standard Izzy STDP 
% FIRST NETWORK CONFIG WITHOUT HOMEOSTASIS: 
% [x] 100 poisson spiking neurons connected to ONE output neuron
% (excitatory regular spiking)
% [x] each poisson neuron had a mean firing rate that ranged from
% [0.2:0.2:20] 
% [x] intial weights S: [aa=.01, bb=.03]  r = aa + (bb-aa).*rand(100,1);
% this is the STDP nearest neighbor implementation: which constrains the
% neuron pairings to only the one before post synaptic  
clear all
close all
clc

Ne=100;       Ni=0;
N = Ne + Ni;
S = zeros(Ne,1);

% parameters of the excitatory RS Izzy neuron 
a= 0.02;
b= 0.2;
c= -65;
d= 8;

% IMPLEMENT 
% LAMBDA is the 

% S contains the weights ordered from, to
% initialize weights according to network1 configuration 
% IMPLEMENT A all-to-one connection 
aa=.01; bb=.03;
for i = 1:N
    for j = 1:1
        S(i,j) = (aa + (bb-aa).*rand(1,1) ) ; %.* 1000;  
    end
end

% STDP lists for LTD and LTP weight changes
LTP(1:Ne, 1:1) = 0;    % ordered from, to
LTD(1:Ne, 1:1) = 0;    % ordered from, to

% STDP parameters
A_plus = 0.0002;
A_minus = 0.000066;
t_plus = 20;
t_minus = 60;

wmin = 0;
wmax = 0.05;

Time=1000;
% there is only 1 v 
v = -65; 
u= b.*v;  % Initial values of u

fr=zeros(Time,1);
spikes=zeros(N,Time);

% keep track of the time steps when the poisson spiked. 
% lambdar=1./ [0.2:0.2:20];

% initialize the first timestep using tau=1/rate 
% for i = 1:N
% poisson_spikes(i,1) = poissrnd( lambdar(i));
% end
% poisson_spikes = zeros(N,Time); 

lambda=[0.2:0.2:20];
for i=1:N
   for t=1:1000-1
      xrand = rand(1);
      % this is generating 'tau' the interspike intervals 
      spikes(i,t+1) =  spikes(i,t) - log(xrand)/lambda(i) ;
   end
end
spikes = ceil(spikes);

for sec = 1:1000  % 1000 simulation seconds 
    firings=0;           % spike timings
    for t=1:1000          % simulation of 1000 ms
       
        I= [4*randn(Ne,1)]; % thalamic input

        fired = find( v>=30 ); % indices of spikes FOR the one IZZY output NEURON

        % Set the input current for all the pre-synaptic neurons that fired
        % very inefficient but clear!!!
        for i = 1:N         % from
                for j=1:1
            I(i) = I(i) + S(i,j) * ( spikes(i,t)==t );  % v(i) >= 30);
                end
        end
        
        if ~isempty(fired)
            firings = [firings; t + 0*fired, fired]; % [timesteo at which it fired, index of neuron that fired]
            disp('fired!')
%             firings = firings+1;
            v(fired) = c(fired);
            u(fired) = u(fired) + d(fired);
            
            % neurons that fired this timestep
            for i=1:size(fired,1)
                
                % if neuron i is excitatory
                % 1)    apply LTP for weights that connect to neuron i 
                %       (i.e., the neurons that fired prior to i)  
                % 2)    apply LTD for weights from neuron i 
                %       (i.e., the neurons that i fired before) 
                % 3)    set the max LTP for  weights that connect from neuron i
                % 4)    set the max LTD for weights that connect to neuron i
                
%                 if fired(i) <= Ne
%                     S(1:Ne, fired(i)) = min(wmax,S(1:Ne, fired(i)) + LTP(:, fired(i))); % LTP to weights that connect to i from all exc 
%                     S(1:Ne, fired(i)) = max(wmin,S(1:Ne, fired(i)) + LTD(fired(i),:)); % LTD to weights that connect to all exc from i
%                     LTP(fired(i),:) = A_plus;   % set max LTP to all exc from i
%                     LTD(:,fired(i)) = A_minus;  % set max LTD to i from all exc
%                 end
                  if fired(i) <= Ne
                    S(1:Ne, fired(i)) = min(wmax,S(1:Ne, fired(i)) + LTP(:, fired(i))); % LTP to weights that connect to i from all exc 
                    S(fired(i),1:Ne) = max(wmin,S(fired(i),1:Ne) + LTD(fired(i),:)); % LTD to weights that connect to all exc from i
                    LTP(fired(i),:) = A_plus;   % set max LTP to all exc from i
                    LTD(:,fired(i)) = A_minus;  % set max LTD to i from all exc
                  end
            end
            
        end;
        v = v + 0.5*(0.04* v.^2 +5*v +140 -u +I);
        v = v + 0.5*(0.04* v.^2 +5*v +140 -u +I);
        u = u + a.*(b.*v-u);
        
        % exponentially decay LTD and LTP based on time constants
        LTP = LTP - LTP/t_plus;
        LTD = LTD - LTD/t_minus;
    end
    
    fr(sec) = firings;
    
%     subplot(1,2,1)
%     if ~isempty(firings)
%         plot(firings(:,1),firings(:,2),'.');
%         title (['Second ', num2str(sec)])
%     end
%     subplot(1,2,2)
%     hist(reshape(S(1:Ne,1:Ne),1,Ne*Ne))
%     pause (0.5)
 
    
end

figure
plot([1:1000], fr) % time, firing rate
figure
plot([1:100],S,'+r')    % synapse ID number, synaptic strength 
    
