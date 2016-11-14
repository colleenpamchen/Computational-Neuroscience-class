%% Midterm_2.m #2 Dopamine modulated STDP solve "distal reward problem"
% (1) Implementation: 
%
% (2) Experimental Design: 
%
% (3) Expected Results: 
%
% Modified DASTDP for Reward PRediction error

% This program employs the use of the programs described below in order to
% use dopamine modulated STDP to make a spiking neural network that will
% account for the probability of recieving rewards given actions and code
% for reward prediction error and uncertainty 

% daspnet.m: Spiking network with DA-modulated STDP. Based on spnet.m
% Created by Eugene M.Izhikevich.                November 3, 2004
%
% 
% Izhikevich E.M. (2007) Solving the Distal Reward Problem through Linkage 
% of STDP and Dopamine Signaling. Cerebral Cortex, 10.1093/cercor/bhl152 

T=3600;
M=50;                 % number of synapses per neuron
D=1;                   % maximal conduction delay 
% excitatory neurons   % inhibitory neurons      % total number 
% Ne=800;                Ni=200;                   N=Ne+Ni;
Ne=400;                Ni=100;                   N=Ne+Ni;
a=[0.02*ones(Ne,1);    0.1*ones(Ni,1)];
d=[   8*ones(Ne,1);    2*ones(Ni,1)];
sm=4;                 % maximal synaptic strength
post=ceil([N*rand(Ne,M);Ne*rand(Ni,M)]); 
post(1:250,1:50)=repmat([251:300],250,1);
s=[ones(Ne,M);-ones(Ni,M)];         % synaptic weights
sd=zeros(N,M);                      % their derivatives
s0 = s(1:50,:);     
s25= s(51:100,:);    
s5= s(101:150,:);
s75=s(151:200,:);
s1=s(201:250,:);
VM=s(251:300,:);
st=15;
%Plotting things
fire=zeros(T,1);
mean1=zeros(T,1);
mean2=zeros(T,1);
mean3=zeros(T,1);
mean4=zeros(T,1);
mean5=zeros(T,1);
followDA = zeros(1000,T); % save DA at each timestep
spikes = zeros(T,1000);
%randomly interleave trials
stimuli=[1,2,3,4,5];%randperm(5);
stim=repmat(stimuli,1, (T/5));
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
firings=[-D 0];                         % spike timings

%---------------
% new stuff related to DA-STDP
T=3600;         % the duration of experiment
DA=0;           % level of dopamine above the baseline
rew=[];

shist=zeros(1000*T,2);
%--------------

for sec=1:T % simulation of 1 day
    %probability distributions for 5 different stimuli
    P=[0,1;.25,.75;.5,.5;.75,.25;1,0];
    %generating reward based on probability
    outcome=mnrnd(1,P(stim(sec),:));
    for t=1:1000                          % simulation of 1 sec
    %saving spike count for seconds and milliseconds
    fire(sec) = fire(sec) + sum(v(251:300)>30);
    spikes(sec,t) = sum((v(251:300)>30));

    %Creating current for all neurons   
    I=13*(rand(N,1)-0.5);               % random thalamic input 
    %show at each time step instead of band
    if t == 250;
        % give current to stimulus neurons such that if there is a 0.0 p of
        % reward give current to neurons 1:50, if p=.25 neurons 51:100 if
        % p=.5 neurons 101-150, if p=.75 neurons 151-200, if p=1 neurons
        % 201-250
     if stim(sec)==1
      I(1:50) = st;    
     elseif stim(sec)==2
      I(51:100) = st;  
     elseif stim(sec)==3
      I(101:150) = st;  
     elseif stim(sec)==4
      I(151:200) = st; 
     elseif stim(sec)==5
      I(201:250) = st;  
      
     end
     %If there is a reward outputted then give dopamine to the system
      if outcome(1,1)==1
%         rew=[rew,sec*1000+2000];
%     end
%     if any(rew==sec*1000+t)
        DA=DA+0.1;
       
    end;
    end
    fired = find(v>=30);                % indices of fired neurons
    v(fired)=-65;  
    u(fired)=u(fired)+d(fired);
    STDP(fired,t+D)=0.5;
    for k=1:length(fired)
      sd(pre{fired(k)})=sd(pre{fired(k)})+STDP(N*t+aux{fired(k)});
    end;
    firings=[firings;t*ones(length(fired),1),fired];
    k=size(firings,1);
    while firings(k,1)>t-D
      del=delays{firings(k,2),t-firings(k,1)+1};
      ind = post(firings(k,2),del);
      I(ind)=I(ind)+s(firings(k,2), del)';
      sd(firings(k,2),del)=sd(firings(k,2),del)-1.5*STDP(ind,t+D)';
      k=k-1;
    end;
    v=v+0.5*((0.04*v+5).*v+140-u+I);    % for numerical 
    v=v+0.5*((0.04*v+5).*v+140-u+I);    % stability time 
    u=u+a.*(0.2*v-u);                   % step is 0.5 ms
    STDP(:,t+D+1)=0.95*STDP(:,t+D);     % .95*STDP tau = 20 ms decay
   
    DA=DA*.95; % decay of dopamine

   
    if (mod(t,10)==0)
        s(1:Ne,:)=max(0,min(sm,s(1:Ne,:)+(0.002+DA)*sd(1:Ne,:)));
        sd=0.9*sd;
    end;
 


% save Dopamine
  followDA(t,sec) = DA;
    
  end;
  %saving mean firing rate
  mean1(sec)=mean(mean((s(1:50,:))));
  mean2(sec)=mean(mean((s(51:100,:))));
  mean3(sec)=mean(mean((s(101:150,:))));
  mean4(sec)=mean(mean((s(151:200,:))));
  mean5(sec)=mean(mean((s(201:250,:))));


% ---- plot -------
 subplot(2,1,1)
  plot(firings(:,1),firings(:,2),'.'); hold on;
  axis([0 1000 0 N]); 
  plot(followDA(:,sec)*N,'r--','linew',5); hold off;
  tit = sprintf('%d',sec);
  title(tit);

%   subplot(2,2,3);
%   plot(fire);
 subplot(2,1,2);
%   hist(s(find(s>0)),sm*(0.01:0.01:1)); % only excitatory synapses
  hold off; 
%   plot(s(1:50,:),repmat(0,50),'r.'); 
%   plot(s(51:100,:),repmat(50,50),'b.');
%   plot(s(101:150,:),repmat(100,50),'m.');
%   plot(s(151:200,:),repmat(150,50),'y.');
%   plot(s(201:250,:),repmat(200,50),'g.');

%hist(spikes(sec,:));

plot(1:T,mean1,'.r'); hold on;
plot(1:T,mean2,'.b')
plot(1:T,mean3,'.m')
plot(1:T,mean4,'.y')
plot(1:T,mean5,'.g')

 % hold off;
% surf(s); shading flat; set(gca,'view',[0   90]);
 drawnow;
  
% ---- end plot ------
  STDP(:,1:D+1)=STDP(:,1001:1001+D);
  ind = find(firings(:,1) > 1001-D);
  firings=[-D 0;firings(ind,1)-1000,firings(ind,2)];
end;
% subplot(5,1,1)
% bar(spikes(501,:),10)
% subplot(5,1,2)
% bar(spikes(502,:),10)
% subplot(5,1,3)
% bar(spikes(503,:),10)
% subplot(5,1,4)
% bar(spikes(504,:),10)
% subplot(5,1,5)
% bar(spikes(505,:),10)
% hold off


