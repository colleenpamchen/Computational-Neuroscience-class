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
alpha = 1; % homeostatic scaling factor? (not given) 
R = 0; % average firing rate of postsynaptic neuron j = OUTPUT neuron 
T = 5000 ; % ms. %T=5 SECONDS. time scale over which the firing rate of the postsunaptic  neuron was averaged 
Rtarget= 35; % 35 Hz target firing rate 
gamma = 50; % homeostasis tuning factor 


Ne=100;       Ni=0;     Nout=1;
N = Ne + Ni;

% parameters of the excitatory RS Izzy neuron 
% a= 0.02;
% b= 0.2;
% c= -65;
% d= 8;
a=[0.02*ones(Ne,1); 0.1*ones(Ni,1)];
b=[0.2*ones(Ne,1);  0.2*ones(Ni,1)];
c=[-65*ones(Ne,1);  -65*ones(Ni,1)];
d=[8*ones(Ne,1);    2*ones(Ni,1)];

% S contains the weights ordered from, to
% initialize weights according to network1 configuration 
% IMPLEMENT A all-to-one connections 
S = zeros(N,1); %% ?????????????? for the weight matrix, is it a
% [poisson neurons 100 x 1 izzy output neuron] ????
% becaues if [N N] then it implies all-to-all connection, which is 
% how izzy had it set up, and later, in the 
aa=.01; bb=.03;
for i = 1:N
    for j = 1:Nout
        S(i,j) = (aa + (bb-aa).*rand(1,1) ) ; %.* 1000;  
    end
end

% STDP lists for LTD and LTP weight changes
LTP(1:Ne, 1:Nout) = 0;    % ordered from, to
LTD(1:Ne, 1:Nout) = 0;    % ordered from, to

% STDP parameters
A_plus = 0.0002;
A_minus = 0.000066;
t_plus = 20;
t_minus = 60;

wmin = 0;
wmax = 0.03;
Time=1000;  % conveninetly, both simulation seconds and timestep ms 

% there is only 1 v 
% v = -65; 
v = -65 * ones( Ne+Ni,1 );
u = b.*v;  % Initial values of u

fr=zeros(Time,1); % FIRING RATE keeps track of how many cycles per second output neuron fired 
spikes=zeros(N,Time); % input POISSON spikes into the system 

% keep track of the time steps when the poisson spiked. 
% lambdar=1./ [0.2:0.2:20];

% initialize the first timestep using tau=1/rate 
% for i = 1:N
% poisson_spikes(i,1) = poissrnd( lambdar(i));
% end
% poisson_spikes = zeros(N,Time); 

lambda=[0.2:0.2:20]; % mean Firing Rate of poisson neurons 
% for i=1:N
%    for t=1:1000-1
%       xrand = rand(1);
%       % this is generating 'tau' the interspike intervals 
%       spikes(i,t+1) =  spikes(i,t) - ( log(xrand)/lambda(i) );
%    end
% end
% spikes = ceil(spikes);

for sec = 1:1000  % 1000 simulation seconds 
    firings=[];           % spike timings
    % this gets reset at every second loop 
    
    
    for t=1:1000          % simulation of 1000 ms
%         I= [4*randn(Ne,1)]; % thalamic input
        I= [40*randn(Ne,1)]; % thalamic input

        fired = find( v>=30 ); % indices of spikes FOR the one IZZY output NEURON
        % original implementation has two dimensions:
        % []
        
        % Set the input current for all the pre-synaptic neurons that fired
        % very inefficient but clear!!!
        for i = 1:N         % from
                for j=1:Nout
            I(i) = I(i) + S(i,j) * (v(i)>= 30);
                end
        end
        
        if ~isempty(fired)
            firings = [firings; t + 0*fired, fired]; % [timestep at which it fired, index of neuron that fired]
            % [timestpe, index of neuron who fired]
%             disp('fired!')
%             fr(sec) = fr+1; % the idea here is to keep track of how
%             many times the izzy output neuron fired to get the avereage output FR
            v(fired) = c(fired);
            u(fired) = u(fired) + d(fired);
            
            % neurons that fired this timestep
            for i=1:size(fired,1)
                
%                 fr(sec) = size(fired,1); % this produced fr around 7 hz, too low! 
                fr(sec)=fr(sec)+1; % the idea here is to keep track of how
                
                R = fr(sec)/100; 
                
                % if neuron i is excitatory
                % 1)    apply LTP for weights that connect to neuron i 
                %       (i.e., the neurons that fired prior to i)  
                % 2)    apply LTD for weights from neuron i 
                %       (i.e., the neurons that i fired before) 
                % 3)    set the max LTP for  weights that connect from neuron i
                % 4)    set the max LTD for weights that connect to neuron i
                

                    S(fired(i) ) = min(wmax,S(fired(i)) + LTP(fired(i))); % LTP to weights that connect to i from all exc 
%                     S(1:Ne, fired(i)) = max(wmin,S(fired(i)) + LTD(fired(i))); % LTD to weights that connect to all exc from i
                    LTP(fired(i)) = A_plus;   % set max LTP to all exc from i
                    LTD(fired(i)) = A_minus;  % set max LTD to i from all exc
                 
            end
            % THIS IS THE WEIGHT CHANGE / UPDATE: APPLY HOMEOSTASIS
            S(fired(i)) = ( alpha * S(fired(i))* (1-(R/Rtarget) ) + ( LTP(fired(i))+LTD(fired(i)) )  ) * (R/ (T*(1+ abs(1-(R/Rtarget)) * gamma) )) ;          
                            
                            
                            
                            
                            
%              fr(sec) = size(fired,1);  % this produced firing rate around
%              6-8hz too low! 
            
        end
        v = v + 0.5*(0.04* v.^2 +5*v +140 -u +I);
        v = v + 0.5*(0.04* v.^2 +5*v +140 -u +I);
        u = u + a.*(b.*v-u);
        
        % exponentially decay LTD and LTP based on time constants
        LTP = LTP - LTP/t_plus;
        LTD = LTD - LTD/t_minus;
    end
    
%     fr(sec) = firings;
    
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
plot((1:Time), fr/100) % time, firing rate
title('HOMEOSTASIS')
figure
plot((1:100),S,'+r')    % synapse ID number, synaptic strength 
title('HOMEOSTASIS')
ylim([-0.01 0.06])    











%% Standard Izzy STDP 
% FIRST NETWORK CONFIG WITHOUT HOMEOSTASIS: 
% [x] 100 poisson spiking neurons connected to ONE output neuron
% (excitatory regular spiking)
% [x] each poisson neuron had a mean firing rate that ranged from
% [0.2:0.2:20] 
% [x] intial weights S: [aa=.01, bb=.03]  r = aa + (bb-aa).*rand(100,1);
% this is the STDP nearest neighbor implementation: which constrains the
% neuron pairings to only the one before post synaptic  
%
% [ ]  FIRST try to make the STDP + Homeostasis learning rule work with
%       non random stimulus.  
% ( ) Without any input POISSON spikes, izzy output neuron fires around 65 hz and weights saturate at .03   (x) PROBLEM: Izzy output neuron isn't firing AT ALL 
% ( ) IMPLEMENT POISSON spikes, and make sure that the regular STDP is
% replicable 
%
%


clear all
close all
clc

Ne=100;       Ni=0;     Nout=1;
N = Ne + Ni;

% parameters of the excitatory RS Izzy neuron 
% a= 0.02;
% b= 0.2;
% c= -65;
% d= 8;
a=[0.02*ones(Ne,1); 0.1*ones(Ni,1)];
b=[0.2*ones(Ne,1);  0.2*ones(Ni,1)];
c=[-65*ones(Ne,1);  -65*ones(Ni,1)];
d=[8*ones(Ne,1);    2*ones(Ni,1)];

% S contains the weights ordered from, to
% initialize weights according to network1 configuration 
% IMPLEMENT A all-to-one connections 
S = zeros(N,1); %% ?????????????? for the weight matrix, is it a
% [poisson neurons 100 x 1 izzy output neuron] ????
% becaues if [N N] then it implies all-to-all connection, which is 
% how izzy had it set up, and later, in the 
aa=.01; bb=.03;
for i = 1:N
    for j = 1:Nout
        S(i,j) = (aa + (bb-aa).*rand(1,1) ) ; %.* 1000;  
    end
end

% STDP lists for LTD and LTP weight changes
LTP(1:Ne, 1:Nout) = 0;    % ordered from, to
LTD(1:Ne, 1:Nout) = 0;    % ordered from, to

% STDP parameters
A_plus = 0.0002;
A_minus = 0.000066;
t_plus = 20;
t_minus = 60;

wmin = 0;
wmax = 0.03;
Time=1000;  % both simulation seconds and timestep ms 

% there is only 1 v 
% v = -65; 
v = -65 * ones( Ne+Ni,1 );
u = b.*v;  % Initial values of u

fr=zeros(Time,1); % FIRING RATE keeps track of how many cycles per second output neuron fired 
spikes=zeros(N,Time); % input POISSON spikes into the system 

% keep track of the time steps when the poisson spiked. 
% lambdar=1./ [0.2:0.2:20];

% initialize the first timestep using tau=1/rate 
% for i = 1:N
% poisson_spikes(i,1) = poissrnd( lambdar(i));
% end
% poisson_spikes = zeros(N,Time); 

lambda=[0.2:0.2:20]; % mean Firing Rate of poisson neurons 
% for i=1:N
%    for t=1:1000-1
%       xrand = rand(1);
%       % this is generating 'tau' the interspike intervals 
%       spikes(i,t+1) =  spikes(i,t) - ( log(xrand)/lambda(i) );
%    end
% end
% spikes = ceil(spikes);

for sec = 1:1000  % 1000 simulation seconds 
    firings=[];           % spike timings
    % this gets reset at every second loop 
    
    
    for t=1:1000          % simulation of 1000 ms
%         I= [4*randn(Ne,1)]; % thalamic input
        I= [40*randn(Ne,1)]; % thalamic input

        fired = find( v>=30 ); % indices of spikes FOR the one IZZY output NEURON
        % original implementation has two dimensions:
        % []
        
        % Set the input current for all the pre-synaptic neurons that fired
        % very inefficient but clear!!!
        for i = 1:N         % from
                for j=1:Nout
            I(i) = I(i) + S(i,j) * (v(i)>= 30);
                end
        end
        
        if ~isempty(fired)
            firings = [firings; t + 0*fired, fired]; % [timestep at which it fired, index of neuron that fired]
            % [timestpe, index of neuron who fired]
%             disp('fired!')
%             fr(sec) = fr+1; % the idea here is to keep track of how
%             many times the izzy output neuron fired to get the avereage output FR
            v(fired) = c(fired);
            u(fired) = u(fired) + d(fired);
            
            % neurons that fired this timestep
            for i=1:size(fired,1)
                
%                 fr(sec) = size(fired,1); % this produced fr around 7 hz, too low! 
                fr(sec)=fr(sec)+1; % the idea here is to keep track of how
                
                % if neuron i is excitatory
                % 1)    apply LTP for weights that connect to neuron i 
                %       (i.e., the neurons that fired prior to i)  
                % 2)    apply LTD for weights from neuron i 
                %       (i.e., the neurons that i fired before) 
                % 3)    set the max LTP for  weights that connect from neuron i
                % 4)    set the max LTD for weights that connect to neuron i
                
%                 if fired(i) <= Ne
                    S(fired(i) ) = min(wmax,S(fired(i)) + LTP(fired(i))); % LTP to weights that connect to i from all exc 
%                     S(1:Ne, fired(i)) = max(wmin,S(fired(i)) + LTD(fired(i))); % LTD to weights that connect to all exc from i
                    LTP(fired(i)) = A_plus;   % set max LTP to all exc from i
                    LTD(fired(i)) = A_minus;  % set max LTD to i from all exc
%                 end
%                   if fired(ii) <= Ne
%                        disp('this should not work!')
%                     S(1:Ne, fired(ii)) = min(wmax,S(1:Ne, fired(ii)) + LTP(:, fired(ii))); % LTP to weights that connect to i from all exc 
%                     S(fired(ii),1:Ne) = max(wmin,S(fired(ii),1:Ne) + LTD(fired(ii),:)); % LTD to weights that connect to all exc from i
%                     LTP(fired(ii),:) = A_plus;   % set max LTP to all exc from i
%                     LTD(:,fired(ii)) = A_minus;  % set max LTD to i from all exc
%                   end
            end
%              fr(sec) = size(fired,1);  % this produced firing rate around
%              6-8hz too low! 
            
        end
        v = v + 0.5*(0.04* v.^2 +5*v +140 -u +I);
        v = v + 0.5*(0.04* v.^2 +5*v +140 -u +I);
        u = u + a.*(b.*v-u);
        
        % exponentially decay LTD and LTP based on time constants
        LTP = LTP - LTP/t_plus;
        LTD = LTD - LTD/t_minus;
    end
    
%     fr(sec) = firings;
    
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
plot([1:Time], fr/100) % time, firing rate
figure
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
