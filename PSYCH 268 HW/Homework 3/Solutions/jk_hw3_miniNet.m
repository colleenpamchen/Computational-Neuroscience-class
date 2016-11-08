% Create a feedforward mean firing rate network.
% Manually tune the weights so that one output neuron responds to vertical 
% and the other responds to horizontal.
function jk_hw3_miniNet

inx = 0;
for i = pi/8:pi/8:pi
    n = runNet(i);
    inx = inx + 1;
    subplot(2,4,inx)
    plot(n)
    axis ([0 25 0 1])
    title (['Angle=', num2str(rad2deg(i))])
    legend('Vertical', 'Horizontal')
    xlabel('timestep')
    ylabel('activity')
end
end


% runNet
%   Input: angle in radians
%   Output: activity of output neurons for T timesteps
function [outN] = runNet (angle)

numInpN = 4;    % number of input neurons
numOutN = 2;    % number of output neurons
VERT = 1;
HORZ = 2;
GAIN = 2;       % gain of the sigmoid function
T = 25;         % number of timesteps

% allocate and initialize neurons
inpN = zeros(1,numInpN);
outN = zeros(T,numOutN);

% allocate and initialize weights
% the input to output weights are set to make vertical directions have
% stronger weights projecting to the vertical output neuron and vice versa
wInpOutExc = zeros(numInpN, numOutN);
pd = pi/numInpN:pi/numInpN:pi;

% used absolute value because it is orientation tuned
%   / - \ |
wInpOutExc(:,VERT) = abs(sin(pd));  
wInpOutExc(:,HORZ) = abs(cos(pd));

% lateral inhibition between output neurons.
wOutOutInh = zeros(numOutN,numOutN);
wOutOutInh(VERT,HORZ) = -2;
wOutOutInh(HORZ,VERT) = -2;

% set initial state of neurons to a random number. just to make the
% dynamics a little more interesting
outN(1,:) = rand(1,numOutN);

for t = 2:T
    
    % set the rates of the input neurons
    for i = 1:numInpN
        inpN(i) = cosTune(pd(i),angle);
    end
    
    % cycle through the output neurons
    for i = 1:numOutN
        
        % calculate the synaptic input based on the presynaptic activity
        % from the previous timestep and the weight from the presynaptic
        % neuron to the postsynaptic neuron
        synIn = 0;
        for j = 1:numInpN
            synIn = synIn + inpN(j) .* wInpOutExc(j,i); % w from j to i
        end
        for j = 1:numOutN
            synIn = synIn + outN(t-1,j) .* wOutOutInh(j,i); % w from j to i
        end
        
        % get the output neuron activity based on the synaptic input
        outN(t,i) = sigmoidN (synIn, GAIN);
    end
end
end

% Cosine Tuning 
% Used absolute value for orientation tuning
%   / - \ |
%   Input:  pd - preferred direction
%           a - angle in radians
%   Output: firing rate based on cosine tuning curve
function r = cosTune(pd, a)
r = abs(cos(a-pd));
end

% Sigmoid Activation
%   Input:  I - synaptic input
%           g - gain
%   Output: firing rate of sigmoid function
function r = sigmoidN (I, g) 
r = 1/(1+exp(-g*I));
end
