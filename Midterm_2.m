%% Midterm_2.m #2 Dopamine modulated STDP solve "distal reward problem"
% (1) Implementation: 
%
% (2) Experimental Design: 
%
% (3) Expected Results: 


M=100;                 % number of synapses per neuron
D=1;                   % maximal conduction delay 
% excitatory neurons   % inhibitory neurons      % total number 
Ne=800;                Ni=200;                   N=Ne+Ni;

% Initialize parameters...
% parapeter ;          % default value;  explanation of parameter
ntrials= 700;         % 400?; number of trials
ntimes = 100;          % 400? or 20ms+10s?; timesteps in each trial
stim = 15;             % current representing stimulus
interval = 20;         % 20; the "coincidence" interval
DA=0;                  % level of dopamine above the baseline
reward_val = .5;           % .5; magnitude fo the increase in DA
punishmag = 0;         % 0;  magnitude of decrease in DA when a2 chosen
DAlr = .8;             % .995; DA decay rate
supnoDA = .002;        % .002; s = s+(supnoDA+DA)*sd
D = 1;                 % 1; max delay (1 -> no delay)
STDP = zeros(N,ntimes+1+D);% initial STDP = 0
ltp = 4;               % .1;  LTP
ltd = -2.5;            % -1.5; LTD (maybe?)
STDPlr=.95;            % .95; STDP decay rate
sddecay=.75;            % .99; sd decay when it is applied
sm=4;                  % max synaptic strength
firings=[-D 0];        % spike timings
ctime = 0.2;           % time constant tau in computing u
v = -65*ones(N,1);     % initial values
u = ctime.*v;          % initial values

a=[0.02*ones(Ne,1);    0.1*ones(Ni,1)]; % STDP paramaters
d=[   8*ones(Ne,1);    2*ones(Ni,1)]; % STDP parameters: extracellular DA 

sm=4;                 % maximal synaptic strength

s=[ones(Ne,M);-ones(Ni,M)];         % synaptic weights
sd=zeros(N,M);                      % their derivatives

sidx =  1:50;       % index for stimulus neurons 
a1idx= 51:100;      % index for action 1 neurons
a2idx=101:150;      % index for action 2 neurons
% index stimulus and action neurons
  sn = s(sidx,:);     % Stimulus neurons
  a1n= s(a1idx,:);    % Action 1 neurons
  a2n= s(a2idx,:);    % Action 2 neurons
  
post=ceil([N*rand(Ne,M);Ne*rand(Ni,M)]); 
post(1:50,1:100) = repmat([a1idx, a2idx], 50, 1);
  
% initialize presynaptic connections with delay structure;
for i=1:N
  if i<=Ne
    for j=1:D
      delays{i,j}=M/D*(j-1)+(1:M/D);
    end;
  else % if inhibitory
    delays{i,1}=1:M;
  end;
  pre{i}=find(post==i & s>0);             % pre excitatory neurons
  aux{i}=N*(D-1-ceil(ceil(pre{i}/N)/(M/D)))+1+mod(pre{i}-1,N);
end;

% initialize things to save;
rewearned= [];  % used to see probability of earning a reward over trials
rew=[];         % Initialize vector of reward delivery times
rewt=[];        % saves reward delivery times w.r.t. trials
avgs = [];      % save average synaptic strength of s->a1 and s->a2
followDA = zeros(ntimes,ntrials); % save DA at each timestep
cumrewvar = nan(ntrials,2); % saves cumulative reward earned


for trial=1:ntrials              % Trial number
  % pick times for stimulus delivery;
  sts = [100];
  startt = 1;  % initialize start and end times for action selection interval
  endt   = 0;  % initialize start and end times for action selection interval
  at     = []; % reset at (used for plotting)
  resp   = []; % reset rewt (used for plotting)
  for t=1:ntimes               % time in ms
    I=13*(rand(N,1)-0.5);      % random current 
    
  % if stimulus presented, give stimulus current to stimulus neurons;
    if any(sts==t)
      I(sidx) = stim;    % give current to stimulus neurons
      startt= t;         % set start t for action selection
      endt = t+interval; % set end time for action selection
      a1s(trial) = 0; a2s(trial) = 0;  % reset action neuron spike counters;
    end
   
  % if in action selection interval, count spikes from a1 and a2 neurons
    if (t>=startt) && (t<=endt)
      a1s(trial) = a1s(trial) + sum(v(a1idx)>30);
      a2s(trial) = a2s(trial) + sum(v(a2idx)>30);
    end
    
    fired = find(v>=30);  % indices of fired neurons
    v(fired)=-65; u(fired)=u(fired)+d(fired); % reset fired neurons
    STDP(fired,t+D)= ltp; % set LTP for fired neurons
    
    % For each neuron that fired, add LTP to the change in weights from its
    %     presynaptic neurons
    for k=1:length(fired)
      sd(pre{fired(k)}) = sd(pre{fired(k)})+STDP(N*t+aux{fired(k)});
    end;
    
  % firings = [time, neurons that fired at this time] (for plotting)
    firings=[firings;t*ones(length(fired),1),fired]; 
    k=size(firings,1); % number of neurons that fired up to this milisecond
    while firings(k,1)>t-D % while kth timestamp is after the longest delay
      del = delays{firings(k,2),t-firings(k,1)+1};
      ind = post(firings(k,2),del); % index of..?
    % add DA modulated current to current I
      I(ind)=I(ind)+s(firings(k,2), del)';
      sd(firings(k,2),del)=sd(firings(k,2),del)+ltd*STDP(ind,t+D)';
      k=k-1; % look at the neurons that fired one timestep earlier
    end;
    
  % update voltage, decay STDP, decay Da concentration;
  % voltage update split for numerical stability
    v=v+0.5*((0.04*v+5).*v+140-u+I); v=v+0.5*((0.04*v+5).*v+140-u+I);
    u=u+a.*(ctime*v-u);               
    STDP(:,t+D+1)=STDPlr*STDP(:,t+D); % decay STDP
    DA=DA*DAlr;                       % decay Da
    
  % update weights every 10 timesteps
    if (mod(t,10)==0)
        s(1:Ne,:)=max(0,min(sm,s(1:Ne,:)+(supnoDA+DA)*sd(1:Ne,:)));
        sd=sddecay*sd;
    end;
    
  % reward delivery;
    if t==endt; % if action selection interval is over
       resp = [resp,t];   % save reponse times on this trial for plotting
    if a1s(trial)>(a2s(trial)); % if a1 chosen
       DA = DA+reward_val; % deliver reward
       rewt = [rewt,t]; % save reward delivery time (for plotting)
       rewearned = [rewearned,1]; % 1 = earned reward.
       rnddelay = 0; %ceil(2000*rand); % random delay to reward
       at   = [at,1]; % save choice for this trial as correct
    else % if wrong action was chosen
       rewearned = [rewearned,0]; % 0 = reward not earned.
       at   = [at,0];   % for plotting
       DA = DA-punishmag; % punish unrewarded firing patterns
    end;
    end;
    
    a1fired=a1s/700; 
    a2fired=a2s/700;
    figure;
    plot(1-a1fired, 'b*') 
    hold on
    plot(1-a2fired,'ro')
    
% save things;
followDA(t,trial) = DA;

  end;
  
  
% % ---- Feedback/Diagnostic plots -------

subplot(2,1,1); hold off;

 sfidx  = ismember(firings(:,2),sidx);  % find stimulus neurons that fired
 a1fidx = ismember(firings(:,2),a1idx); % find a1 neurons that fired
 a2fidx = ismember(firings(:,2),a2idx); % find a2 neurons that fired

 plot(firings(:,1),firings(:,2),'k.'); hold on; % plot all firings
 plot(firings(a1fidx,1),firings(a1fidx,2),'r.');% plot a1  firings
 plot(firings(a2fidx,1),firings(a2fidx,2),'b.');% plot a2  firings
 plot(firings(sfidx,1),firings(sfidx,2),'g.');  % plot s   firings

 plot([0,1000],[max(sidx),max(sidx)],'g-'); hold on;

 plot([0,1000],[max(a1idx),max(a1idx)],'r-');
 plot([0,1000],[max(a2idx),max(a2idx)],'b-'); % add plots for rewad as a red x or green o at the reward time and over the action that produces reward
 plot(firings(:,1),firings(:,2),'k.');

ylabel('Index of Fired Neuron')
xlabel('Time (ms)')

 plot(resp(at==1),at(at==1),'go') % correct response
 plot(resp(at==0),at(at==0),'rs') % incorrect response

 plot(followDA(:,trial)*N,'k--')

tit = sprintf('trial %d',trial);
title(tit)
axis([0 ntimes 0 N]); % x is time, y is spikes

 subplot(2,2,4)
 hold off; 
 plot(s(sidx,a1idx-50),'.r'); hold on; % corrrct
 plot(s(sidx,a2idx-50),'.b')  % incorrect
 plot([1,50],mean(mean(s(sidx,a1idx-50))).*ones(2,1),'r--') % correct
 plot([1,50],mean(mean(s(sidx,a2idx-50))).*ones(2,1),'b--') % incorrect
 title('Synaptic Strength')
 ylabel('Synaptic Strength')
 xlabel('Stimulus Neuron Index')
 ylim([0 sm])
 subplot(2,2,4)
 hold off; 
 surf(s); shading flat; set(gca,'view',[0,90]);
 colorbar;
 zlim([-1,4])
 title('Synaptic Strength')
 ylabel('Presynaptic Neuron')
 xlabel('Postsynaptic Neuron')

cumrewvar(trial,1) = sum(rewearned)./ntrials;   % compute cumulative gain
cumrewvar(trial,2) = sum(rewearned==0)./ntrials;% compute cumulative loss


subplot(2,2,3); hold on;
plot(1:ntrials,cumrewvar(:,1),'b-') % plot rewards gained
plot(1:ntrials,cumrewvar(:,2),'r-') % plot rewards lost

ylim([0,1])
xlim([0,ntrials])

title('Model Performance')
ylabel('Proportion of Total Possible Rewards')
xlabel('Trial')

drawnow;
% ---- end plotting ------

  STDP(:,1:D+1)=STDP(:,(ntimes+1):(ntimes+1+D));
  ind = find(firings(:,1) > (ntimes+1)-D);
  firings=[-D 0;firings(ind,1)-ntimes,firings(ind,2)];

end;


%%
% STDP = zeros(N,1001+D);
% v = -65*ones(N,1);                      % initial values
% u = 0.2.*v;                             % initial values
% firings=[-D 0];                         % spike timings
% 
% %---------------
% % new stuff related to DA-STDP
% T=3600;         % the duration of experiment
% DA=0;           % level of dopamine above the baseline
% rew=[];
% 
% n1=1;           % presynaptic neuron
% syn=1;          % the synapse number to the postsynaptic neuron
% n2 = post(n1,syn) % postsynaptic neuron
% s(n1,syn)=0;    % start with 0 value
% 
% interval = 20;  % the coincidence interval for n1 and n2
% n1f=-100;       % the last spike of n1
% n2f=[];         % the last spike of n2
% shist=zeros(1000*T,2);
% %--------------
% 
% for sec=1:T                             % simulation of 1 day
%   for t=1:1000                          % simulation of 1 sec
%     I=13*(rand(N,1)-0.5);               % random thalamic input 
%     fired = find(v>=30);                % indices of fired neurons
%     v(fired)=-65;  
%     u(fired)=u(fired)+d(fired);
%     STDP(fired,t+D)=0.1;
%     for k=1:length(fired)
%       sd(pre{fired(k)})=sd(pre{fired(k)})+STDP(N*t+aux{fired(k)});
%     end;
%     firings=[firings;t*ones(length(fired),1),fired];
%     k=size(firings,1);
%     while firings(k,1)>t-D
%       del=delays{firings(k,2),t-firings(k,1)+1};
%       ind = post(firings(k,2),del);
%       I(ind)=I(ind)+s(firings(k,2), del)';
%       sd(firings(k,2),del)=sd(firings(k,2),del)-1.5*STDP(ind,t+D)';
%       k=k-1;
%     end;
%     v=v+0.5*((0.04*v+5).*v+140-u+I);    % for numerical 
%     v=v+0.5*((0.04*v+5).*v+140-u+I);    % stability time 
%     u=u+a.*(0.2*v-u);                   % step is 0.5 ms
%     STDP(:,t+D+1)=0.95*STDP(:,t+D);     % tau = 20 ms
%    
%     DA=DA*0.995;
%     if (mod(t,10)==0)
%         s(1:Ne,:)=max(0,min(sm,s(1:Ne,:)+(0.002+DA)*sd(1:Ne,:)));
%         sd=0.99*sd;
%     end;
%     if any(fired==n1)
%         n1f=[n1f,sec*1000+t];
%     end
%     if any(fired==n2)
%         n2f=[n2f,sec*1000+t];
%         if (sec*1000+t-n1f(end)<interval) & (n2f(end)>n1f(end))
%             rew=[rew,sec*1000+t+1000+ceil(2000*rand)];
%         end;
%     end     
%     if any(rew==sec*1000+t)
%         DA=DA+0.5;
%     end;
%     shist(sec*1000+t,:)=[s(n1,syn),sd(n1,syn)];
% 
%   end;
% % ---- plot -------
%   subplot(2,1,1)
%   plot(firings(:,1),firings(:,2),'.');
%   axis([0 1000 0 N]); 
%   
%   subplot(2,2,3);
%   plot(0.001*(1:(sec*1000+t)),shist(1:sec*1000+t,:), 0.001*rew,0*rew,'rx');
%  
%   subplot(2,2,4);
%   hist(s(find(s>0)),sm*(0.01:0.01:1)); % only excitatory synapses
%   hold on; 
%   
%   plot(s(n1,syn),0,'r.'); hold off;
%   drawnow;
% % ---- end plot ------
%   STDP(:,1:D+1)=STDP(:,1001:1001+D);
%   ind = find(firings(:,1) > 1001-D);
%   firings=[-D 0;firings(ind,1)-1000,firings(ind,2)];
% end;
% 
% 
