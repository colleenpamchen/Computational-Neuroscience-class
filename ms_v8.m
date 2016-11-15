% Notes from paper;
% 1000 neurons
% 50 stimulus, 50 a1, 50 a2.
% spikes measured for 20ms after stimulus presentation
% 400 trials, seperated by 10s
% so.. stimulus, 20ms interval counting spikes, 10s, stimulus, ...

% Instrumental STDP
%   "stimuli" neurons, "action 1" neurons, "action 2" neurons.
%    Later, add response chains, ratio schedules (default case is a FR1
%    schedule, change to FRn or variable ratio) and interval schedules, add
%    more complex relationships between actions and rewards. Add a stimuli
%    for reward features (so instead of just increasing DA, have a set of
%    neurons represent the rewarding outcome).

% Debugging;
    % Try reducing connectivity (down from 100)
    % Try some kind of normalization on the weights to enhance competition
    % between the neurons.
    % change reward size (DA), LTD (1.5), or decay on STDP
    % debugging to make sure reward is being delivered when it should be..
    % First, try this; Try changing connections to be random.. 
    
    % Diagnostic plots for DA concentration, STDP heat map...
    % For heat maps, plot a dot on the y-axis of action neurons (red for
    % a1, blue for a2, black for s)? Maybe plot a1,a2,s clusters
    % seperately with a 'o' in the cells that correspond to a connection
    % with an action or stimulus neuron?

% Resources;
%   Sequence Learner
%   Polychronous groups; Or Bucci, Chou, Krichmar (2014). 
% Send-fire chain (less complicated..). (neurons can represent stimuli or
% acitions..)
%   Zheng & Trisch 2014 Robust development of synfire chains from multiple
%   plasticity mechanisms.

% So..
% 1) initialize neurons 80e/20i
% 2) Pick a random group of excitatory neurons to represent a stimuli
% 3) Pick two random groups of excitatory neurons to represent responses
% 4) Stimuli excited responses, if response A activity > response B, a
%      reward increases DA.
% ^ Sounds good!

% Steps;
% 1) Choose 50 stimulus (S), 50 action 1, and 50 action 2 neurons.
% 1b) Make sure S is connected to a1 and a2 neurons!
% 2) Inject current into S neurons
% 3) Count spikes from a1 and a2 clusters for 20ms after stimuli.
% 4) for 400 trials, reward if a1 spikes strictly more than a2
% 5) for 400 trials, reward if a2 spikes strictly more than a1
% 6) 

% later; 
%   try a reversal learning network? 
%   So a pool for S1->A1, a pool for S1->A2, a pool for S1->A1, and S2->A2?


%% Begin script;
clear;

% Diagnostic plotting & feedback (1 = on, 0 = off); 
%   note: plotting significantly slows things down;
dotback    = 1; % feedback dots (cheap way to know code is still running)
figure;         % initialize a figure for plotting 
alab       = 0; % label axes (slows things down a little)
plotsyn    = 1; % plot synaptic strengths after each trial
plotcumrew = 0; % plot and compute 'online' cumulative rewards gained/lost
plotraster = 1; % plot raster
   impresp = 0; %      impose response (go = correct, ro = incorrect)
   imprew  = 0; %      impose reward delivery (b* = reward)
   impDA   = 1; %      impose DA level across time

% Initialize options...
% note options start with "opt" and should be 0 or 1;
% option ;      % default value;  explanation of option
optswt0 = 1;    % 1; 1 = initialize stimulus->neuron weights as 0?
optstoa = 1;    % 1; 1 = force s neurons to connect to action neurons
optrndsa= 0;    % 0; 1 = choose random s and a neurons

% Initialize neurons and network...
% excitatory   % inhibitory   % total       % connection per neuron
  Ne=800;        Ni=200;        N=Ne+Ni;       M = 200; 
% s = synaptic weight matrix, rows = "from neuron", col = "to neuron"
s=[ones(Ne,M);-ones(Ni,M)];  % synaptic weights; (e:1, i:-1)
sd=zeros(N,M);               % initialize weight change matrix

% Initialize parameters...
% parapeter ;          % default value;  explanation of parameter
ntrials= 700;         % 400?; number of trials
ntimes = 400;          % 400? or 20ms+10s?; timesteps in each trial
stim = 15;             % current representing stimulus
interval = 20;         % 20; the "coincidence" interval
DA=0;                  % level of dopamine above the baseline
rewmag = .5;           % .5; magnitude fo the increase in DA
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
d=[   8*ones(Ne,1);    2*ones(Ni,1)];   % STDP parameters

% pick 50 stimulus, 50 a1, and 50 a2 neurons...
if optrndsa == 1; % randomly
sidx = randsample(1:Ne,50,false);                       % Stimulus neurons
a1idx= randsample(setdiff(1:Ne,sidx),50,false);         % Action 1 neurons
a2idx= randsample(setdiff(1:Ne,[sidx a1idx]),50,false); % Action 2 neurons
else % not randomly
sidx =  1:50;       % index for stimulus neurons 
a1idx= 51:100;      % index for action 1 neurons
a2idx=101:150;      % index for action 2 neurons
% index stimulus and action neurons
end
  sn = s(sidx,:);     % Stimulus neurons
  a1n= s(a1idx,:);    % Action 1 neurons
  a2n= s(a2idx,:);    % Action 2 neurons
  
% randomize post synaptic connections,
% post is # of neurons by # of connections,
% post = random excitatory connections to any 100 neurons, random
%        inhibitory connections to any of the excitatory neurons
post=ceil([N*rand(Ne,M);Ne*rand(Ni,M)]); 
if optstoa ==1;% if force stimulus neurons to connect to a1 and a2 neurons;
post(1:50,1:100)=repmat([a1idx,a2idx],50,1);
end

% initialize presynaptic connections with delay structure;
for i=1:N
  if i<=Ne % if excitatory
    for j=1:D
      %delays{neuron,delay}=(1:connections/delays+connections*delay/delays)
      delays{i,j}=M/D*(j-1)+(1:M/D);
    end;
  else % if inhibitory
      % delays{neuron,1} = 1:connections
    delays{i,1}=1:M;
  end;
  % presynaptic neurons for neuron i are the excitatory (s>0) ones that are
  %  connected to i (post == i)?
  pre{i}=find(post==i&s>0);
  aux{i}=N*(D-1-ceil(ceil(pre{i}/N)/(M/D)))+1+mod(pre{i}-1,N);
end;

if optswt0 == 1; % if initialize s-> neuron weights as 0;
s(1:50,:)=zeros(50,M); % 0; Stimulus -> neuron weights    
end

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
       DA=DA+rewmag; % deliver reward
       rewt = [rewt,t]; % save reward delivery time (for plotting)
       rewearned = [rewearned,1]; % 1 = earned reward.
       rnddelay = 0;%ceil(2000*rand); % random delay to reward
       at   = [at,1]; % save choice for this trial as correct
    else % if wrong action was chosen
       rewearned = [rewearned,0]; % 0 = reward not earned.
       at   = [at,0];   % for plotting
       DA = DA-punishmag; % punish unrewarded firing patterns
    end;
    end;
    
% save things;
followDA(t,trial) = DA;

  end;
  
  
% % ---- Feedback/Diagnostic plots -------
if dotback == 1 % feedback dots
if mod(trial,10)==0; fprintf('.'); end
if mod(trial,100)==0; fprintf('%d/%d, \n',trial,ntrials); end
end

if plotraster == 1
subplot(2,1,1); hold off;
if optrndsa == 1; % if s,a1,a2 neurons are randomly assigned;
 sfidx  = ismember(firings(:,2),sidx);  % find stimulus neurons that fired
 a1fidx = ismember(firings(:,2),a1idx); % find a1 neurons that fired
 a2fidx = ismember(firings(:,2),a2idx); % find a2 neurons that fired
 plot(firings(:,1),firings(:,2),'k.'); hold on; % plot all firings
 plot(firings(a1fidx,1),firings(a1fidx,2),'r.');% plot a1  firings
 plot(firings(a2fidx,1),firings(a2fidx,2),'b.');% plot a2  firings
 plot(firings(sfidx,1),firings(sfidx,2),'g.');  % plot s   firings
else              % if s,a1,a2 indexed in order
 plot([0,1000],[max(sidx),max(sidx)],'g-'); hold on;
 if optstoa == 1; %    and if s forced to connect to a's
 plot([0,1000],[max(a1idx),max(a1idx)],'r-');
 plot([0,1000],[max(a2idx),max(a2idx)],'b-'); % add plots for rewad as a red x or green o at the reward time and over the action that produces reward
 end
 plot(firings(:,1),firings(:,2),'k.');
end

if alab ==1 % label axes
ylabel('Index of Fired Neuron')
xlabel('Time (ms)')
end

if impresp == 1; % if impose response,
 plot(resp(at==1),at(at==1),'go') % correct response
 plot(resp(at==0),at(at==0),'rs') % incorrect response
end

if imprew == 1; % if impose reward delivery,
 if ~isempty(rewt)
 plot(rewt,at,'k*')
 rewt=[];
 end
end

if impDA == 1 ; % if impose DA concentration;
 plot(followDA(:,trial)*N,'k--')
end

tit = sprintf('trial %d',trial);
title(tit)
axis([0 ntimes 0 N]); % x is time, y is spikes
end % end raster plot

% plot synaptic strength
if plotsyn == 1; % if plot synapses
 if (optstoa == 1) && (optrndsa==0); % if s->a'a forced, and a's indexed
 subplot(2,2,4)
 hold off; 
 plot(s(sidx,a1idx-50),'.r'); hold on; % corrrct
 plot(s(sidx,a2idx-50),'.b')  % incorrect
 plot([1,50],mean(mean(s(sidx,a1idx-50))).*ones(2,1),'r--') % correct
 plot([1,50],mean(mean(s(sidx,a2idx-50))).*ones(2,1),'b--') % incorrect
 if alab ==1 % label axes
 title('Synaptic Strength')
 ylabel('Synaptic Strength')
 xlabel('Stimulus Neuron Index')
 end
 ylim([0 sm])
 else % if s->neurons random
 subplot(2,2,4)
 hold off; 
 surf(s); shading flat; set(gca,'view',[0,90]);
 colorbar;
 zlim([-1,4])
 if alab ==1 % label axes
 title('Synaptic Strength')
 ylabel('Presynaptic Neuron')
 xlabel('Postsynaptic Neuron')
 end
 end
end

if plotcumrew == 1; % if plot cumulative reward
cumrewvar(trial,1) = sum(rewearned)./ntrials;   % compute cumulative gain
cumrewvar(trial,2) = sum(rewearned==0)./ntrials;% compute cumulative loss
subplot(2,2,3); hold on;
plot(1:ntrials,cumrewvar(:,1),'b-') % plot rewards gained
plot(1:ntrials,cumrewvar(:,2),'r-') % plot rewards lost
ylim([0,1])
xlim([0,ntrials])
if alab ==1 % label axes
title('Model Performance')
ylabel('Proportion of Total Possible Rewards')
xlabel('Trial')
end
end

%  subplot(2,2,3);
% plot(0.001*(1:(trial*1000+t)),shist(1:trial*1000+t,:), 0.001*rew,0*rew,'rx');
%  subplot(2,2,4);
%  hist(s(find(s>0)),sm*(0.01:0.01:1)); % only excitatory synapses
% %   hold on; plot(s(n1,syn),0,'r.'); hold off;

% plot connections
% subplot(1,1,1)
% surf(s); shading flat; set(gca,'view',[0   90]); 
% ylabel('From neuron'); xlabel('To neuron');

drawnow;
% ---- end plotting ------


  STDP(:,1:D+1)=STDP(:,(ntimes+1):(ntimes+1+D));
  ind = find(firings(:,1) > (ntimes+1)-D);
  firings=[-D 0;firings(ind,1)-ntimes,firings(ind,2)];

end;

%% Model Evaluation Plotting; 

figure;
csum = cumsum(rewearned);
plot(cumsum(rewearned)./sum(rewearned)); hold on;
h = legend('Prop.'); 
set(h,'loc','best');
ylabel('Porportion of total rewards earned by this trial')
xlabel('trial')
title('Rate of total rewards accumulated')

figure;
plot(cumsum(rewearned)./ntrials); hold on;
plot([0 ntrials],[0 1]);
plot([0 ntrials],[0 .5]);
plot([0 ntrials],[sum(rewearned) sum(rewearned)]./ntrials,'k--');
h = legend('Prop.','Always correct', 'guessing', 'total rewards earned'); 
set(h,'loc','best');
ylabel('Porportion of total possible rewards accumulated')
xlabel('trial')
title('Rate of total possible rewards accumulated')


