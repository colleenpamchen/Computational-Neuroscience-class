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
% 
% 
%% Standard Izzy STDP 
% FIRST NETWORK CONFIG WITHOUT HOMEOSTASIS: 
% [] 100 poisson spiking neurons connected to an excitatory regular spiking
% output neron
% [] each poisson neuron had a mean firing rate that ranged from
% [0.2:0.2:20] ???????how to implement/set this
% [x] intial weights S: [aa=.01, bb=.03]  r = aa + (bb-aa).*rand(100,1);
% 

Ne=100;                Ni=0;
N = Ne + Ni;

a=[0.02*ones(Ne,1); 0.1*ones(Ni,1)];
b=[0.2*ones(Ne,1);  0.2*ones(Ni,1)];
c=[-65*ones(Ne,1);  -65*ones(Ni,1)];
d=[8*ones(Ne,1);    2*ones(Ni,1)];

aa=.01; bb=.03;
% S contains the weights ordered from, to
% initialize weights according to network1 configuration 
for i = 1:N
    for j = 1:N
        if i <= Ne
            S(i,j) = aa + (bb-aa).*rand(1,1);
        else
            S(i,j) = -1.0;
        end
    end
end

% STDP lists for LTD and LTP weight changes
LTP(1:Ne, 1:Ne) = 0;    % ordered from, to
LTD(1:Ne, 1:Ne) = 0;    % ordered from, to

% STDP parameters
A_plus = 0.0002;
A_minus = 0.000066;
t_plus = 20;
t_minus = 60;

wmin = 0;
wmax = 0.05;

v=-65*ones(Ne+Ni,1);  % Initial values of v
u=b.*v;               % Initial values of u

for sec = 1:600  % ??????? 1000 simulation seconds 
    firings=[];           % spike timings
    for t=1:1000          % simulation of 1000 ms
        I=[4*randn(Ne,1);2*randn(Ni,1)]; % thalamic input

        fired=find(v>=30); % indices of spikes

        % Set the input current for all the pre-synaptic neurons that fired
        % very inefficient but clear!!!
        for i = 1:N         % from
            for j = 1:N     % to
                I(j) = I(j) + S(i,j) * (v(i) >= 30);
            end
        end
        
        if ~isempty(fired)
            firings=[firings; t+0*fired, fired];
            v(fired)=c(fired);
            u(fired)=u(fired)+d(fired);
            
            % neurons that fired this timestep
            for i=1:size(fired,1)
                
                % if neuron i is excitatory
                % 1)    apply LTP for weights that connect to neuron i 
                %       (i.e., the neurons that fired prior to i)  
                % 2)    apply LTD for weights from neuron i 
                %       (i.e., the neurons that i fired before) 
                % 3)    set the max LTP for  weights that connect from neuron i
                % 4)    set the max LTD for weights that connect to neuron i
                
                if fired(i) <= Ne
                    S(1:Ne, fired(i)) = min(wmax,S(1:Ne, fired(i)) + LTP(:, fired(i))); % LTP to weights that connect to i from all exc 
                    S(fired(i),1:Ne) = max(wmin,S(fired(i),1:Ne) + LTD(fired(i),:)); % LTD to weights that connect to all exc from i
                    LTP(fired(i),:) = A_plus;   % set max LTP to all exc from i
                    LTD(:,fired(i)) = A_minus;  % set max LTD to i from all exc
                end
            end
            
        end;
        v=v+0.5*(0.04*v.^2+5*v+140-u+I);
        v=v+0.5*(0.04*v.^2+5*v+140-u+I);
        u=u+a.*(b.*v-u);
        
        % exponentially decay LTD and LTP based on time constants
        LTP = LTP - LTP/t_plus;
        LTD = LTD - LTD/t_minus;
    end
    
    subplot(1,2,1)
    if ~isempty(firings)
        plot(firings(:,1),firings(:,2),'.');
        title (['Second ', num2str(sec)])
    end
    subplot(1,2,2)
    hist(reshape(S(1:Ne,1:Ne),1,Ne*Ne))
    pause (0.5)
end
