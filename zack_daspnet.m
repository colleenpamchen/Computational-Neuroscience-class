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
rng('default')
M=100;                 % number of synapses per neuron
D=1;                   % maximal conduction delay
% excitatory neurons   % inhibitory neurons      % total number
Ne=800;                Ni=200;                   N=Ne+Ni;
a=[0.02*ones(Ne,1);    0.1*ones(Ni,1)];
d=[   8*ones(Ne,1);    2*ones(Ni,1)];
sm=4;                 % maximal synaptic strength
 
%generates a list of connections for each neuron
post=ceil([N*rand(Ne,M);Ne*rand(Ni,M)]);
s=[ones(Ne,M);-ones(Ni,M)];         % synaptic weights
sd=zeros(N,M);                      % their derivatives
%for all neurons
for i=1:N
    if i<=Ne %if this is an excitatory neuron
        for j=1:D %D = 1
            delays{i,j}=M/D*(j-1)+(1:M/D);
        end;
    else  %if this is an inhibitory neuron
        delays{i,1}=1:M;
    end;
    %delays will always be 1:100 for this particular case
     
    % list of indices of presynaptic connections for each neuron
    pre{i}=find(post==i&s>0);             
    % this extremely convoluted calculation finds the inverse list of
    % presynaptic neurons for each neuron.  For example, if neuron 1 has
    % neuron 592 as a presynaptic neuron, you will see -408 appear in its
    % list in aux{}.  All you do is add 1000 and you will get the correct
    % neuron index
    aux{i}=N*(D-1-ceil(ceil(pre{i}/N)/(M/D)))+1+mod(pre{i}-1,N);
end;
%these are STDP values for all the neurons for 1 second plus 1ms plus D
%second delay, which in this case is 1ms
STDP = zeros(N,1001+D);
v = -65*ones(N,1);                      % initial values
u = 0.2.*v;                             % initial values
firings=[-D 0];                         % spike timings
 
%---------------
% new stuff related to DA-STDP
T=3600;         % the duration of experiment
DA=0;           % level of dopamine above the baseline
rew=[];
 
n1=1;           % presynaptic neuron
syn=1;          % the synapse number to the postsynaptic neuron
n2=post(n1,syn) % postsynaptic neuron
s(n1,syn)=0;    % start with 0 value
 
interval = 20;  % the coincidence interval for n1 and n2
n1f=-100;       % the last spike of n1
n2f=[];         % the last spike of n2
shist=zeros(1000*T,2);
%--------------
 
for sec=1:T   
    sec% simulation of 1 day
    for t=1:1000                          % simulation of 1 sec
        I=13*(rand(N,1)-0.5);               % random thalamic input
        fired = find(v>=30);                % indices of fired neurons
        v(fired)=-65;
        u(fired)=u(fired)+d(fired);
        % for every neuron that fired, set the STDP value of the next time
        % step (t + D where D=1) to 0.1
        STDP(fired,t+D)=0.1;
        %for every neuron that fired
        for k=1:length(fired)
            %for each synapse it has,
            %the derivative of the synaptic weight for the presynaptic
            %connections to the fired neuron increase by the STDP value at
            %the current timestep
            sd(pre{fired(k)})=sd(pre{fired(k)})+STDP(N*t+aux{fired(k)});
        end;
        %firings is a record of:
        %   timestep, nueron that fired
        %where there can be multiple of one timestep such as:
        %   1   14
        %   1   15
        %   1   16
        %meaning neurons 14, 15, and 16 all fired at timestep 1.
        firings=[firings;t*ones(length(fired),1),fired];
        k=size(firings,1);
        %do for every firing that happened at a timestep greater than this
        %one minus D (which equals 1)
        %in other words, do this for every neuron that fired this timestep
        while firings(k,1)>t-D
            %del is the delay(neuron that fired, delay list number, which
            %is always 1, since firings(k,1) is always t, so this is really
            %t-t+1
            %if D were greater than one this could conceviably be a
            %different list of delays, but in this case there is only one
            %list of delays, which is always 1 to 100
            del=delays{firings(k,2),t-firings(k,1)+1};
            %get the indices of the post synaptic neurons for the one that
            %fired
            ind = post(firings(k,2),del);
            %the current for those post synaptic neurons increases by s,
            %the synaptic weight
            I(ind)=I(ind)+s(firings(k,2), del)';
            %the derivative of these weights (synapse weights for current 
            %neuron that fired this timestep for all its 100 synapses)
            %changes by -1.5 times the STDP value of these neurons, 
            %which is either .1 if it fired this step or 
            %0.1 *.95^number of timesteps since it fired
            sd(firings(k,2),del)=sd(firings(k,2),del)-1.5*STDP(ind,t+D)';
            k=k-1;
        end;
        %voltage stuff
        v=v+0.5*((0.04*v+5).*v+140-u+I);    % for numerical
        v=v+0.5*((0.04*v+5).*v+140-u+I);    % stability time
        u=u+a.*(0.2*v-u);                   % step is 0.5 ms
        STDP(:,t+D+1)=0.95*STDP(:,t+D);     % tau = 20 ms
        %value of STDP decays each timestep, so next timestep's STDP values
        %are 95% of this timestep's STDP values
         
        %level of extracellular dopamine decays each timestep
        DA=DA*0.995;
        %every 10 timesteps
        if (mod(t,10)==0)
            %sets synaptic strength, limited to 0 <= s <= sm(which is 4)
            %0.002 + DA is baseline DA amount (2nM) + amount of DA
            %so it is DA amount * rate of change in synaptic strength
            %Remember, this is for all synapses
            s(1:Ne,:)=max(0,min(sm,s(1:Ne,:)+(0.002+DA)*sd(1:Ne,:)));
            %decay rate of change of synaptic strength
            %eligibility trace in paper
            %play with this decay value (increase for longer reward distances)      
            sd=0.99*sd;
        end;
 
    %Look at short term plasticity, pgs 184-185 5.37
    %example on class website for CARL roomba, stp()
 
 
        %for tracking the target neuron
        if any(fired==n1)
            n1f=[n1f,sec*1000+t];
        end
        if any(fired==n2)
            n2f=[n2f,sec*1000+t];
            %if neuron 2 fired within 20ms after neuron 1
            %set a reward to happen randomly 1 to 3 seconds in the future
            if (sec*1000+t-n1f(end)<interval) & (n2f(end)>n1f(end))
                rew=[rew,sec*1000+t+1000+ceil(2000*rand)];
            end;
        end
        %if there is a reward this timestep
        %increase the extracellular DA by 0.5
        if any(rew==sec*1000+t)
            DA=DA+0.5;
        end;
        shist(sec*1000+t,:)=[s(n1,syn),sd(n1,syn)];
         
    end;
    % ---- plot -------
    subplot(2,1,1)
    plot(firings(:,1),firings(:,2),'.');
    axis([0 1000 0 N]);
    subplot(2,2,3);
    plot(0.001*(1:(sec*1000+t)),shist(1:sec*1000+t,:), 0.001*rew,0*rew,'rx');
    subplot(2,2,4);
    hist(s(find(s>0)),sm*(0.01:0.01:1)); % only excitatory synapses
    hold on; plot(s(n1,syn),0,'r.'); hold off;
    drawnow;
    % ---- end plot ------
    STDP(:,1:D+1)=STDP(:,1001:1001+D);
    ind = find(firings(:,1) > 1001-D);
    firings=[-D 0;firings(ind,1)-1000,firings(ind,2)];
end;
