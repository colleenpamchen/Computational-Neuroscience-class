% %% Midterm_2.m #2 Dopamine modulated STDP solve "distal reward problem"
% (1) Implementation: 
%
% (2) Experimental Design: 
%
% (3) Expected Results: 
% daspnet.m: Spiking network with DA-modulated STDP. Based on spnet.m
% Created by Eugene M.Izhikevich.                November 3, 2004
%
% This program reproduces the experiment in Fig.1 in
% Izhikevich E.M. (2007) Solving the Distal Reward Problem through Linkage 
% of STDP and Dopamine Signaling. Cerebral Cortex, 10.1093/cercor/bhl152 
%
% n1 - the presynaptic neuron. syn is the synapse to be reinforced.
% Plot: top - spike raster. Bottom left - synaptic strength (blue), the
% eligibility trace (green), and the rewards (red x). Bottom right - the
% distribution of synaptic weights with the chosen synapse marked by red dot.
clear all
close all
clc 

M=100;                 % number of synapses per neuron
D=1;                   % maximal conduction delay 
% excitatory neurons   % inhibitory neurons      % total number 
Ne=800;                Ni=200;                    N=Ne+Ni;
% Ne=400;                Ni=100;                   N=Ne+Ni;
a=[0.02*ones(Ne,1);    0.1*ones(Ni,1)];
d=[   8*ones(Ne,1);    2*ones(Ni,1)];
sm=4;                 % maximal synaptic strength
I = zeros(N,1); 

randidx = randperm(Ne,150);
S = randidx(1:50);   
A = randidx(51:100);
B = randidx(101:150);


post=ceil([N*rand(Ne,M);Ne*rand(Ni,M)]); 
s=[ones(Ne,M);-ones(Ni,M)];         % synaptic weights
sd=zeros(N,M);                      % their derivatives
for i=1:N
  if i<=Ne
    for j=1:D
      delays{i,j}=M/D*(j-1)+(1:M/D);
    end;
  else
    delays{i,1}=1:M;
  end;
  pre{i}=find(post==i&s>0);             % pre excitatory neurons
  aux{i}=N*(D-1-ceil(ceil(pre{i}/N)/(M/D)))+1+mod(pre{i}-1,N);
end;
STDP = zeros(N,1001+D);
v = -65*ones(N,1);                      % initial values
u = 0.2.*v;                             % initial values
% a list of timestamps and what neuron fired. intialized at D=-1 because
% nothing fired then 
firings=[-D 0];                         % spike timings

%---------------
% new stuff related to DA-STDP
T=800;         % number of trials
DA=0;           % level of dopamine above the baseline


% for neurons A and B, find the presynaptic neurons that are in group S
interval = 20;  % the coincidence interval for n1 and n2
% shist=zeros(10000*T,2);
%--------------

probA=0;
probB=0;

for trial=1:T     % simulation of 1 trial
    rew=0;
    
  for t=1:10000    % simulation of 1 msec
    % stimulate S neurons for 1ms at start
    if t==1
        I(S) = 30;
    else
        I(S) = 0;
    end
    %--------don't touch--------------
    %find neurons that fired, reset their v and u, set STDP
    fired = find(v>=30);                % indices of fired neurons
    v(fired)=-65;  
    u(fired)=u(fired)+d(fired);
    STDP(fired,t+D)=0.1;
    for k=1:length(fired)
      sd(pre{fired(k)})=sd(pre{fired(k)})+STDP(N*t+aux{fired(k)});
    end;
    firings=[firings;t*ones(length(fired),1),fired];
    %----------------------------------
    %------------don't touch-----------
    % record of firing history as [timestamp, neuron]
    % can have multiple of timestamp
    k=size(firings,1);
    while firings(k,1)>t-D
      del=delays{firings(k,2),t-firings(k,1)+1};
      ind = post(firings(k,2),del);
      I(ind)=I(ind)+s(firings(k,2), del)';
      sd(firings(k,2),del)=sd(firings(k,2),del)-1.5*STDP(ind,t+D)';
      k=k-1;
    end;
    %----------------------------------
    
    % TODO
    if t==20
          % using firings,
          % count the number of neurons in group A that fired |A|
          numA = sum( ismember(firings(:,2),A ) );
          % count the number of neurons of group B that fired |B|
          numB = sum( ismember(firings(:,2),B ) );
          % compare & reward
         if numB > numA && trial >= 401 
            % set rew = timestamp random in the future up to 1s, inversely
            % proportional to |A|/|B| or |B|/|A| depending on which trial set
            rew = t + (numB/numA)*1000 ; 
         end
    end
        
    
    %----------don't touch---------
    v=v+0.5*((0.04*v+5).*v+140-u+I);    % for numerical 
    v=v+0.5*((0.04*v+5).*v+140-u+I);    % stability time 
    u=u+a.*(0.2*v-u);                   % step is 0.5 ms
    STDP(:,t+D+1)=0.95*STDP(:,t+D);     % tau = 20 ms
   
    DA=DA*0.995;
    if (mod(t,10)==0)
        s(1:Ne,:)=max(0,min(sm,s(1:Ne,:)+(0.002+DA)*sd(1:Ne,:)));
        sd=0.99*sd;
    end;
    %------------------------------
    
    if rew == t % if reward is now, give reward 
        DA=DA+0.5;
    end;
    
    %TODO figure out history for plotting
%     shist(trial*10000+t,:)=[s(n1,syn),sd(n1,syn)];

% calculate probabilities of response for A or response for B : 
% foreach trial
probA(trial)=  numA/(numA+numB); disp(probA)
probB(trial)=  numB/(numA+numB); disp(probB)
% plot()



  end; % end of t loop 
  
% ---- plot -------
  subplot(2,1,1)
  plot(firings(:,1),firings(:,2),'.');
  axis([0 1000 0 N]); 
  subplot(2,2,3);
%   plot(0.001*(1:(trial*1000+t)),shist(1:trial*1000+t,:), 0.001*rew,0*rew,'rx');
  subplot(2,2,4);
  hist(s(find(s>0)),sm*(0.01:0.01:1)); % only excitatory synapses
  hold on; plot(s(n1,syn),0,'r.'); hold off;
  drawnow;
% ---- end plot ------
  STDP(:,1:D+1)=STDP(:,1001:1001+D);
  ind = find(firings(:,1) > 1001-D);
  firings=[-D 0;firings(ind,1)-1000,firings(ind,2)];
  
end;  % END OF time loop 

