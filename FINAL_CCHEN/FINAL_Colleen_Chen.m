%% FINAL EXAM: Colleen Chen 
%  
% [ ]  FIRST try to make the STDP + Homeostasis learning rule work with
%       non random stimulus.  
% ( ) alpha- homeostatic scaling factor??  
% ( ) R average firing rate of the postsynpatic neuron j

clear all
close all
clc

% ADD NEW homeostatic variables: 
alpha = 1; % homeostatic scaling factor? (not given) 
R = 36; % average firing rate of postsynaptic neuron j = OUTPUT neuron 
T = 5; %5000 ; % ms. %T=5 SECONDS. time scale over which the firing rate of the postsunaptic  neuron was averaged 
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

% STDP lists for LTD and LTP weight changes
LTP1(1:Ne, 1:Nout) = 0;    % ordered from, to
LTD1(1:Ne, 1:Nout) = 0;    % ordered from, to

% STDP parameters
A_plus = 0.0002;
A_minus = 0.000066;
t_plus = 20;
t_minus = 60;

Time  =      1000;  % conveninetly, both simulation seconds and timestep ms 

wmin = 0;
wmax = 0.041;

v1 = -65 ; % there is only 1 v
u1 = b*v1;  % Initial values of u

% S contains the weights ordered from, to
% initialize weights according to network1 configuration 
% IMPLEMENT A all-to-one connections 
S1 = zeros(N,1); 
aa=.01; bb=.03;
for i = 1:N
    for j = 1:Nout
        S1(i,j) = (aa + (bb-aa).*rand(1,1) ) ; 
    end
end

fr1 = zeros(Time,1); % FIRING RATE keeps track of how many cycles per second output neuron fired 
I1=zeros(N,1); 

% [ ] keep track of the time steps when the poisson spiked. 
% [ ] LAMBDA would be the value of pixels, all 28*28 values. W= (64) izzy neuron
% would have 784 weights [784, 64]
% H= firing rate of 64 neurons  [64, 1] 
% [784,64] * [64,1] = [784, 1] image 


count1 = zeros(N,1);
lambda=[0.2:0.2:20]'; % mean Firing Rate of poisson neurons per SEC   

spikes1=zeros(N,1);
xrand = rand(N,1); 
spikes1(:,2) =  spikes1(:,1) - ( log(xrand)./lambda ) ; % this is generating 'tau' the interspike intervals

for sec = 1:Time  % 1000 simulation seconds 
%     disp(sec);
    dsdt = zeros(N,1);
    if sec<2
        R=36;
    else
    R = fr1(sec);
    end
    
    K = R / ( T.*( 1+ abs(1-R/Rtarget)*gamma ) ) ;

    for t=1:1000        % simulation of 1000 ms   
        % IF POISSON SPIKED  
         Ifired1 = spikes1(:,2) <= t/1000+sec ; % logical array that allows you to index 
        
    if ~isempty(Ifired1) 
            spikes1(Ifired1,1) = spikes1(Ifired1,2);
            xrand = rand(sum(Ifired1),1); 
            spikes1(Ifired1,2) = spikes1(Ifired1,1) - (log(xrand)./lambda(Ifired1)) ;     
    
            pidx1 = find(Ifired1); %         index of all poisson i's that spiked:

            count1(Ifired1) = count1(Ifired1)+1;   % total number of times each poisson neuron fired 
            numFired1 = sum(Ifired1);              % number of poisson neurons that fired at this time step
            % update voltage
            I1(Ifired1) = I1(Ifired1)+ S1(Ifired1); 

                if v1 >=30 
                   fr1(sec)=fr1(sec)+1;
                   v1 = c; % reset V 
                   u1 = u1 + d;

                   % LTP to those that fired, LTD to all else 
%                    dsdt = dsdt+ ( alpha .* S1 .* (1 - R/Rtarget) + 1.*(LTP1 + LTD1) ).* K ;
                     dsdt = dsdt+ K.* ( alpha .* S1 .* (1 - R/Rtarget) ) + K.*(LTP1 + LTD1)  ;
%                    dsdt = dsdt+ (alpha .* S1 .* (1 - R/Rtarget) + 1.*(LTD1)).* K;
                   LTP1(Ifired1) = A_plus;
%                    LTD1(Ifired1) = A_minus;

                end % if V fired...     
                dsdt = dsdt+ K.*(alpha .* S1 .* (1 - R/Rtarget) ) + K.*(LTD1) ;
                LTD1(Ifired1) = A_minus;
        end % end of Ifired1 
        v1 = v1 + 0.5*( 0.04* v1^2 +5*v1 +140 -u1 + sum(I1(Ifired1)) ) ;
        v1 = v1 + 0.5*(0.04* v1^2 +5*v1 +140 -u1 + sum(I1(Ifired1)) ) ;
        u1 = u1 + a.*(b.*v1-u1);

        LTP1 = LTP1 - LTP1/t_plus;
        LTD1 = LTD1 - LTD1/t_minus;
    end % end of t ms loop 

        S1 = S1 + dsdt;  
%         S1 = max(wmin, S1 ) ;
        S1 = max(wmax, S1 ) ;
  
   
end % End of TIME second loop 
poisson_firing_rate_given1 = count1/Time;

% subplot(2,2,1)
% plot((1:Time), fr1/28) % time, firing rate
% title('HOMEOSTASIS')
% % figure
% subplot(2,2,2)
% plot((1:100),S1,'+r')    % synapse ID number, synaptic strength 
% title('HOMEOSTASIS')
% ylim([-0.001 0.04]) 


% .......................................................

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
% spikes = zeros(N,Time); % input POISSON spikes into the system 
Ifired=[];
I=zeros(N,1);  

% INITIALIZE poisson, simulate 5 sec of poison spikes:
spikes = zeros(N,2); % input POISSON spikes into the system 
count = zeros(N,1);

lambda=[0.2:0.2:20]'; % mean Firing Rate of poisson neurons per SEC   
xrand=rand(N,1); 
spikes(:,2) =  spikes(:,1) - ( log(xrand)./lambda ) ; % this is generating 'tau' the interspike intervals

for sec = 1:Time  % 1000 simulation seconds 
    disp(sec);
    firings = []; 
    for t=1:1000        % simulation of 1000 ms   
        % IF POISSON SPIKED    
%         t/1000+sec 
        Ifired = spikes(:,2) <= t/1000+sec ; % logical array that allows you to index 
        
        if ~isempty(Ifired) 
        firings = [firings; t+0*Ifired,Ifired];
        spikes(Ifired,1) = spikes(Ifired,2);
        xrand = rand(sum(Ifired),1); 
        spikes(Ifired,2) = spikes(Ifired,1) - (log(xrand)./lambda(Ifired)) ;     
        pidx = find(Ifired); % index of all poisson i's that spiked      
        count(Ifired) = count(Ifired)+1; 
        numFired = sum(Ifired);
%             str = [num2str( numFired), ' poisson fired at ', num2str(t) ];
%             disp(str);           
           
            % update voltage
            I(Ifired) = I(Ifired)+ S(Ifired); 
       
            if v >=30 
               fr(sec)=fr(sec)+1;
               v = c; % reset V 
               u = u + d;
               S(Ifired) = min(wmax, S(Ifired) + LTP(Ifired)); % LTP to weights that connect to all exc from i
                % S = min(wmax, S + LTP); % LTP to weights that connect to all exc from i
               LTP(Ifired) = A_plus;  % set max LTD to i from all exc
                %             S(~Ifired) = max(wmin, S(~Ifired)+ LTD(~Ifired)); % LTP to weights that connect to i from all exc  
                % % S(Ifired) = max(wmin, S(Ifired)+ LTD(Ifired));
                %             LTD(Ifired) = A_minus;   
            end
            S(Ifired) = max(wmin, S(Ifired)+ LTD(Ifired)); % LTP to weights that connect to i from all exc  
            % S(Ifired) = max(wmin, S(Ifired)+ LTD(Ifired));
            LTD(Ifired) = A_minus;

end % end of Ifired1 
            v = v + 0.5*(0.04* v^2 +5*v +140 -u + sum(I(Ifired)) ) ;
            v = v + 0.5*(0.04* v^2 +5*v +140 -u + sum(I(Ifired)) ) ;
            u = u + a.*(b.*v-u);
            
            LTP = LTP - LTP/t_plus;
            LTD = LTD - LTD/t_minus;
end % end of t ms loop 
  
        S = max(wmin, S);    
   
end % End of TIME second loop 
poisson_firing_rate_given = count/Time;

% subplot(2,2,3)
% plot([1:Time], fr ) % time, firing rate
% subplot(2,2,4)
% plot([1:100],S,'+r')    % synapse ID number, synaptic strength 
% % ylim([-0.01 0.06]) 

subplot(1,2,2)
plot((1:Time), fr1 ,'r') % time, firing rate
title('Firing rate of output neuron with and without Homeostasis enabled')
hold on
plot([1:Time], fr ,'b') % time, firing rate
hold off
%...
subplot(1,2,1)
plot([1:100],S,'+b')    % synapse ID number, synaptic strength 
title('Synaptic strength of synapses with and without Homeostasis enabled')
hold on
% ylim([-0.01 0.06])  
plot((1:100),S1,'+r')    % synapse ID number, synaptic strength 
hold off
% ylim([-0.001 0.04]) 


%% Standard Izzy STDP 
% FIRST NETWORK CONFIG WITHOUT HOMEOSTASIS: 
% IMPLEMENT A all-to-one connections: 

% S = zeros(N,1); % [poisson neurons 100 x 1 izzy output neuron] 
% aa=.01; bb=.03;
% for i = 1:N
%         S(i) = (aa + (bb-aa).*rand(1,1) ) ; % Synaptic weights 100 x 1 outputNeuron 
% end
% 
% % STDP lists for LTD and LTP weight changes
% LTP(1:Ne, 1:Nout) = 0;    % ordered from, to
% LTD(1:Ne, 1:Nout) = 0;    % ordered from, to
% 
% % there is only 1 v 
% v = -65; 
% u = b*v;  % Initial values of u
% 
% fr = zeros(Time,1); % FIRING RATE keeps track of how many cycles per second output neuron fired 
% spikes = zeros(N,Time); % input POISSON spikes into the system 
% fired=[];
% I=zeros(N,1); 
% % nspikes = zeros(N,1);
% 
% lambda=[0.2:0.2:20]; % mean Firing Rate of poisson neurons 
% for i=1:N
%    for t=1:1000 % I guess this is looping through the Time trials from 1:1000 SECONDS
%       xrand = rand(1);
%       % this is generating 'tau' the interspike intervals 
%       spikes(i,t+1) =  spikes(i,t) - ( log(xrand)/lambda(i) ) ; % generate 1000 spikes per neuron in SECONDs  
%    end
% end
% spikes = ceil(spikes); % these are the spike times in ms.  
% 
% for sec = 1:Time  % 1000 simulation seconds 
%     disp(sec)
%     vfired=[]; 
%     
%     for t=1:1000 %1000          % simulation of 1000 ms
%        
%         for i = 1:N  % this loops through the poisson input NEURONS 
%           
%             if v >=30
%                 fr(sec)=fr(sec) + 1;
%                 v = c;
%                 u = u + d;   
%                 vfired = [vfired; t]; % v neuron fired at t 
% %                  disp('vfired!')
%                 S(i) = min(wmax, S(i) + LTP(i)); % LTP to weights that connect to all exc from i
%                 LTP(i) = A_minus;  % set max LTD to i from all exc
%             end 
%             
%                 Ifired = find(spikes(i,1:Time+1)==t); % the t indices of ith poisson that fired 
%                 if ~isempty(Ifired) 
%                         I(i) = I(i) + S(i) ;%* S(Ifired);  %*length(fired) ; % If at this timestep, neuron i spiked, then add input and synaptic weight
%                        
%                         S(i) = max(wmin, S(i)+ LTD(i)); % LTP to weights that connect to i from all exc 
%                         LTD(i) = A_plus;   % set max LTP to all exc from i          
%                 end
%              v = v + 0.5* (0.04 * v^2 + 5*v + 140 -u + I(i) );
%              v = v + 0.5*(0.04 * v^2 + 5*v + 140 -u + I(i) );       
%              u = u + a*(b*v -u ); 
%         end % poissonNeuron i loop
%         LTP = LTP - LTP/t_plus;
%         LTD = LTD - LTD/t_minus;
%     end % end of t ms loop    
%     
% end % end of Time s loop 
% 
% % figure
% subplot(2,2,3)
% plot([1:Time], fr./5 ) % time, firing rate
% % figure
% subplot(2,2,4)
% plot([1:100],S,'+r')    % synapse ID number, synaptic strength 
% ylim([-0.01 0.06])    


% 
% %%
% % Neuron network setup:
% %   100 poisson spiking neurons
% %   Excitatory Izhekevich regular spiking (RS)
% %   Inhibitory fast spiking (FS)
% %   delay: 1ms 
% %   initial weights: randn()/uniform distribution [.01 .03]
% %   STDP paramters: 
% %       apre, apost, taupre, taupost 
% % 
% % 
% % synaptic weights
% % firing rates
% % target homeostatsis firing rate: 35hz 
% % 
% % Homeostasis & STDP Learning rule:   
% % 
% % weight update
% % 
% % 
% % 
% % 
% % poisson_spikes = zeros(N,1000); 
% % poisson_spikes(1,1)=0;
% % for n=1:100
% %    for t=1:1000-1
% %        xrand = rand(1);
% %       poisson_spikes(n,t+1) =  poisson_spikes(n,t) - log(xrand)/lambda(i) ;
% % %       fired = find(poisson_spikes>xrand); 
% %    end
% % end
% % 
% %  
% % % 
% % % 
% % % % 
% % for i = 1:N
% % poisson_spikes(i,:) = poissrnd( lambda(i),1,1000 );
% % end
% % 
% % % 
% % 
