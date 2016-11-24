% Psych268A - Fall 2016
% Midterm Exam
%   Demontrate divisive normalization in a small network
function jk_midterm_divisive_normalization

% Preferred Directions
global pd;
pd = [ -pi/2 -pi/4 -pi/8 0 pi/8 pi/4 pi/2 ];

% Baseline Firing Rate
global r0;
r0 = 36;

% Calculate the cosine tuning curves for visualization purposes
for i = 1:size(pd,2)
    inx = 0;
    for j = -pi:0.1:pi
        inx = inx + 1;
        r(inx) = cosineTuning (j, r0, 60, pd(i));
    end
    plot(r)
    hold on
end

% Give a range of maximum firing rates to simulate contrast. Low rMax
% represents low contrast and high rMax represents high contrast
rMax = [40, 60, 80, 100 ];

figure;

% Replicate Figure 3e from Carandini and Heeger
%   Top row has one stimulus
stimulus = 0;
for i = 1:size(rMax,2)
    fr = runDivisiveNormalization(rMax(i), stimulus, 1);
    subplot(2,size(rMax,2),i)
    plot(fr)
    axis([0 size(fr,2)+1 0 1])
    title(['Div Norm - Rmax = ', num2str(rMax(i))])
end

%   Bottom row has two stimuli
stimulus = [pi/2; -pi/2];
for i = 1:size(rMax,2)
    fr = runDivisiveNormalization(rMax(i), stimulus, 1);
    subplot(2,size(rMax,2),i+size(rMax,2))
    plot(fr)
    axis([0 size(fr,2)+1 0 1])
    title(['Div Norm - Rmax = ', num2str(rMax(i))])
end

figure;

% Re-run the model without normalization
stimulus = 0;
for i = 1:size(rMax,2)
    fr = runDivisiveNormalization(rMax(i), stimulus, 0);
    subplot(2,size(rMax,2),i)
    plot(fr)
    axis([0 size(fr,2)+1 0 1])
    title(['No Norm - Rmax = ', num2str(rMax(i))])
end

stimulus = [pi/2; -pi/2];
for i = 1:size(rMax,2)
    fr = runDivisiveNormalization(rMax(i), stimulus, 0);
    subplot(2,size(rMax,2),i+size(rMax,2))
    plot(fr)
    axis([0 size(fr,2)+1 0 1])
    title(['No Norm - Rmax = ', num2str(rMax(i))])
end

end

function fr = runDivisiveNormalization(rMax, stimulus, divisiveNormalization)

global r0
global pd;
numInput = size(pd,2);
numOutput = numInput;
numInhibitory = 1;
nInput = zeros(1,numInput);
nOutput = zeros(1,numOutput);
nInhibitory = zeros(1,numInhibitory);
rHighContrast = 100;

gain = 0.1;

% If there is not divisive normalization set the weights lower to get good
% activity levels. Using the same weights maxed out activity for all
% neurons.
if divisiveNormalization
    W_INP_OUT = 1;
else
    W_INP_OUT = 0.10;
end
wInpOut = zeros(numInput,numOutput);
for i = 1:numInput
    wInpOut(i,i) = W_INP_OUT;
end

% There is a many to one projection from the input neurons to the inhibitory
% neuron. So, the weight is relatively weak.
W_INP_INH = 0.002;
wInpInh = zeros(numInput,numInhibitory);
for i = 1:numInput
    for j = 1:numInhibitory
        wInpInh(i,j) = W_INP_INH;
    end
end

% There is a one to many projection from the inhibitory neuron to the output neurons.
% So the weight is negative and strong.
W_INH_OUT = -125;
wInhOut = zeros(numInhibitory,numOutput);
for i = 1:numInput
    for j = 1:numInhibitory
        wInhOut(i,j) = W_INH_OUT;
    end
end


T = 10;
for ts = 1:T
    
    % Make a copy of the previous activity. Use this when calculating
    % activity
    nInhibitoryPrev = nInhibitory;
    nInputPrev = nInput;
    
    % Get the firing rate, based on a cosine tuning curve, for each input neuron.
    for i = 1:numInput
        
        % If there are 2 stimuli, the second stimulus is a high contrast
        % input.
        if size(stimulus,1) == 2
            nInput(i) = cosineTuning(stimulus(1), r0, rMax, pd(i)) + cosineTuning(stimulus(2), r0, rHighContrast, pd(i));
        else
            nInput(i) = cosineTuning(stimulus(1), r0, rMax, pd(i));
        end
    end
    
    % Get the inhibitory neuron activity. The synaptic input is from the 
    % input neuron activity times the weights. Use a sigmoid function to get the
    % activity. If there is no divisive normalization, set the inhibitory
    % activity to zero.
    for i = 1:numInhibitory
        if (divisiveNormalization)
            nInhibitory(i) = sigmoid( sum(nInputPrev .* wInpInh(i,:)) , gain);
        else
            nInhibitory(i) = 0;
        end
    end
    
    % Get the output neuron activity. The synaptic input comes from the
    % input layer and the inhibitory neuron. Run this through a sigmoid
    % function.
    for i = 1:numOutput
        nOutput(i) = sigmoid( sum( nInput .* wInpOut(i,:)) + sum(nInhibitoryPrev .* wInhOut(i,:)) , gain);
    end
end

fr = nOutput;

disp(['Rmax = ', num2str(rMax), ' Divisive Normalization = ', num2str(divisiveNormalization)])
fprintf(1,'   Input: ');
fprintf(1,'%3.3f ', nInput); fprintf('\n')
fprintf(1,'   Inhibitory: ');
fprintf(1,'%3.3f ', nInhibitory); fprintf('\n')
fprintf(1,'   Output: ');
fprintf(1,'%3.3f ', nOutput); fprintf('\n')
fprintf('\n')

end

function s = sigmoid(synIn, gain)

s = 1/(1 + exp(-gain*synIn));

end

function f = cosineTuning (input, r0, rMax, sMax)

f = max(0, r0 + (rMax - r0) * cos(input - sMax));

end