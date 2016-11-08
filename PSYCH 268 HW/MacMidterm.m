function MacMidterm
close all; 

plots = 0; % make 0 to turn off the plot of S matrix by trial

problem1;
% 1) implimentation
%    I initialize weights for the two conditions, then provide each with
%    the same input but block input for the ttx condition. See Problem1
%    function for a more detailed description.
% 2) design
%    Weights update according to a simple Hebbian rule.
% 3) expected results
%    I expected results that qualitatively matched the figures in the
%    guiding paper.

problem2(plots);
% 1) implimentation
%    I started with the issy script and made many modifications to the
%    connection structure, then tested it and did not get output that
%    looked like the article.. I then used numerous diagnostic plots to
%    inform changes to the structure and parameters until I got good
%    output.
%    See the Problem2 function below for more details.
% 2) design
%    100 CA3 neurons connect to one CA1 neuron. 
% 3) expected results
%    I expectd the results to qualitatively match the figures in the
%    guiding paper.


end

% Problem 1
function problem1
% conditions; sutured (s) & ttx inactivated (ttx)
% Results; Due to long term depression (LTD), greater reduction in response
%          was seen in the sutured eye (which had spontanious activity)
%          than the ttx condition (which had no activation).
% Measures; OD = proportion of a neuron's response that is driven by
%                stimulation of the right vs left eye. 
%              = (right - left) / (right + left)

clear;
% initialize
e = .01; % learning rate
ts = 0:100; % time steps

% weights = [right eye, left eye]
% 100x2; 100 neurons by 2 input neurons
ws = 2.4*randn(100,2);  % weights for ttx
wttx = ws;              % weights for sutured

% Train the weights;
for t = ts;
% for sutured
inps = rand(2,1); % random input
dws = e.*ws*inps*inps'; % simple hebb rule..
%dw = outp * (inp-(sumv*inp/n))'; % hebb + normalization?
% c = (inp-mean(inp))*(inp-mean(inp))'; ws = w*c; % cov rule?
ws = ws + dws; % update weights for sighted group

% for ttx
inpt = inps.*[1;0]; % block input to right eye.
dwt = e.*wttx*inpt*inpt'; % simple hebb rule..
%dw = outp * (inpb-(sumv*inpb/n))'; % hebb + normalization?
% c = (inp-mean(inp))*(inp-mean(inp))'; wb = w*c; % cov rule?
wttx = wttx + dwt; % updateweights for blind group
end

% outputs in each eye [right,left] for sighted and blind
outs = ws; 
outt = wttx;

% Occular dominance
ODs = max(min((outs*[1;-1])./(outs*[1;1]),1),-1);
ODt = max(min((outt*[1;-1])./(outt*[1;1]),1),-1);

% Plotting
figure(1); hold on;
% histograms
hx = subplot(2,2,1);
histogram(sort(ODt),5,'normalization','prob', ...
          'BinWidth',.4,'BinLimits',[-1,1],...
          'facecolor',[0,0,0]);
title('Sutured'); ylim([0 .5])
set(hx,'xtick',-.8:.4:.8); set(hx,'xticklabel',1:5)
ylabel('Neurons (%)'); xlabel('DE < ---------- > OD')
hx = subplot(2,2,2);
histogram(sort(ODs),5,'normalization','prob', ...
          'BinWidth',.4,'BinLimits',[-1,1],...
          'facecolor',[.7,.7,.7]);
title('TTX'); ylim([0 .5])
set(hx,'xtick',-.8:.4:.8); set(hx,'xticklabel',1:5)
ylabel('Neurons (%)'); xlabel('DE < ---------- > OD')

% Cumulative distribution of OD
subplot(2,1,2); hold on;
plot(sort(ODt),1:100,'-.','color',[0,0,0])
plot(sort(ODs),1:100,'-.','color',[.6,.6,.6])
h = legend('Sutured','TTX'); set(h,'location','best')
title('Cumulative Distribution of Neuron OD')
xlabel('Occular dominance, -1 = deprived eye, 1 = open eye')
ylabel('Number of neurons')

% left and right activation
% figure;
% subplot(3,2,5)
% plot(sort(outs))
% title('TTX')
% h = legend('left activation','right activation'); set(h,'location','best')
% ylabel('Output'); xlabel('Sorted Neuron Index')
% subplot(3,2,6)
% plot(sort(outt))
% title('Sutured')
% h = legend('left activation','right activation'); set(h,'location','best')
% ylabel('Output'); xlabel('Sorted Neuron Index')


end





% problem 2
function problem2(plots)

% Notes:
% 1. CA1 neuron is the 101st neuron in S matrix.

plotting = plots; % plotting on if plots == 1


% initialize CA3s
nc3 = 100; % 100 c3 neurons
c3max = 40; % max input for CA3-->CA1 connection

% initialize CA1 weights
CA1w = normpdf(1:100,50,10);
CA1w = CA1w./max(CA1w)*4;

% initialize S matrix; S(from neuron, to neuron)
S = zeros(nc3+1); S(:,end)=[CA1w 0]; 

% number of neurons;
ns = length(S);

% initialize LTD and LTP weight changes
LTP(1:ns, 1:ns) = 0; % LTP(from, to)
LTD(1:ns, 1:ns) = 0; % LTD(from, to)

% initialize Izzy & STDP parameters and variables
a=1.02*ones(ns,1); %.02*ones(nc3+1,1); originally
b=0.2*ones(ns,1);
c=-65*ones(ns,1);
d=8*ones(ns,1);
maxLTP  = .1; 
maxLTD = -.105;
Pdecay =  10;        % 20; originally
Ddecay =  20;
wmin = 0;
wmax = 5.0;
v=-65*ones(nc3+1,1); % Initial values of v (potential)
vthresh = 30;        % threshold for spiking
u=b.*v;              % Initial values of u (recovery variable)

% initialize loop variables
nlaps= 17; nlocs = 100;
laps = 1:nlaps;   % rat runs 17 laps.
locs = 1:nlocs;  % rat experiences 100 locations each lap.

% initialize saving variables;
CAf = zeros(nlaps,nlocs); % CAf(lap,loc) saves CA1 fired (1) or not (0)

for lap = laps % for each lap...
for loc=locs % Step through locations
  % set location dependent CA3 inputs;
    I = zeros(ns,1);                  % reset inputs
    I(1:nc3)= normpdf(loc,1:nc3,3.9); % gaussian tuning for CA3 neurons
    I = I./max(I)*c3max;              % rescale to CA3 max rate
    fired=find(v>=vthresh);         % indices of neurons that fired
    CAf(lap,loc) = any(fired==101); % 1 if CA1 fired, 0 otherwise
    CA3f = setdiff(fired,101);      % indecies of CA3 neurons that fired
    I(end) =  sum(S(CA3f,end)); % CA1 input = sum of weights on fired CA3s
  %  hold off; plot(I);pause(.01); % diagnostic plot for I
    
  % if >=1 neuron has fired, Update connection matrix (S) 
    if ~isempty(fired)     % if >=1 neuron has fired
        v(fired)=c(fired); % reset fired neurons to resting potential
        u(fired)=u(fired)+d(fired); 
        
      % for each neuron that fired at this location 
        for ii=1:size(fired,1)
      % if CA1 neuron fired
            if fired(ii) == 101;
          % Potentiate all connections to CA1
            S(:, 101) = min(wmax,S(:, 101) + LTP(:, 101)); 
          % set max LTD on all connections to CA1
            LTD(:,101) = maxLTD;  
            end
      % if CA3 neurons fired
          if fired(ii) <= nc3
          % Depotentiate connection from fired CA3 to CA1
            S(fired(ii),101) = max(wmin,S(fired(ii),101) + LTD(fired(ii),101)); 
          % set LTP on all connections from fired CA3s
            LTP(fired(ii),:) = maxLTP;
          end
         end 
    end;

% plot S if user specified plots == 1;
if plotting == 1;
figure(2)
hold off; plot(loc,I(loc),'k*'); hold on; 
surf(S); shading interp;
set(gca,'view',[-24.3000   41.2000]);
tit = sprintf('Lap: %d/17',lap);
title(tit);
zlim([0 5])
drawnow;
end

% plot S for only CA3 --> CA1 neuron weights
% plot(S(:,end)); hold on; plot(loc,2,'k*'); hold off; drawnow;

   % compute voltage
   v=v+0.5*(0.04*v.^2+5*v+140-u+I);
   v=v+0.5*(0.04*v.^2+5*v+140-u+I);
  % hold off; plot(v); pause(.01); % diagnostic plot for v 
  
  % compute recovery var
   u =u+a.*(b.*v-u);
  % plot(u); pause(.01) % diagnostic plots for u
      
  % exponentially decay LTD and LTP based on time constants
    LTP = LTP*(1-1/Pdecay);
    LTD = LTD*(1-1/Ddecay);
      

end
% save skewness after each lap;
  CA1Inskew(lap) = skewness(S(:,end));

  
end

% plotting;
% Fig a.
% 1. Initial strength symmetric, magnitude proportional to red lines in (a).
% 2. grey lines -> 0 weight
% 3. red curve = net input before learning 
% 4. CA1 also receives oscillatory inhibitory input
% 5. Synaptic matrix becomes asymmetric (blue lines and blue curve) after experience.
figure;
subplot(2,2,1); hold on;
plot(CA1w,'red')
plot(S(1:nc3,end),'blue')
h = legend('Inithal CA1 input','Experienced CA1 input'); set(h,'location','best')
ylabel('Strength')
xlabel('Location')
title('Fig. a')

% Fig b.
subplot(2,2,3); hold on;
plot(conv(CAf(1,:),normpdf(-1:.1:1)),'r-')
plot(conv(CAf(end,:),normpdf(-1:.1:1)),'b-')
legend('First lap','Last lap')
xlabel('Location')
ylabel('Firing rate (arbitrary units)')
title('Fig b.')
% 1. red = first lap CA1 activity
% 2. blue = CA1 activity after experience
 
% Fig d.
subplot(2,2,2); hold on;
plot(CA1Inskew);
legend('input');
ylabel('skewness'); xlabel('lap')
title('Fig d.')
% 1. green = skewness of input of CA3-> CA1 synaptic matrix
% 2. blue = skewness of CA1 output

% Fig e.
subplot(2,2,4); hold on;
% find first and last spikes;
firstf = []; lastf = [];
for ii = 1:size(CAf,1)
firstf = [firstf find(CAf(ii,:),1,'first')];
lastf = [lastf find(CAf(ii,:),1,'last')];
end
% 1. Blue: location where first CA1 spike was observed
plot(firstf,'b')
% 2. Green: location where last CA1 spike observed
plot(lastf,'g')
h = legend('first spike','last spike'); set(h,'location','best')
xlabel('lap')
ylabel('location')
title('Fig e.')



end





