%% midterm #1 Divisive Normalization   author: Colleen Chen
% (1) Implementation: 
% The input neuron activities are calcuated by fitting the gaussian tuning curve
% through the input stimulus with each of the neuron's preferred orientation. Then, for the input layer neurons,
% activities (fr) are multiplied by each of the two sets of synaptic weights: (1) excitatory, 
% (2) inhibitory. The input to output layer is of an excitatory connection,
% so the dot product of the positive weights and the firing rates from
% input neurons, fed through a sigmoid activation function, gives us the
% excitatory output activities. The inhibitory layer, where the divisive normalization occurs,
% actually sums up the dot product of the input neurons and inhibitory weights, and then another
% layer of inhibition with dot product of negative inhibitory weights are
% fed through a sigmoid activation function for inhibitory output
% activities. Sum up the inhibitory and exchitatory output activities. We
% have, for each neuron, a summary activity. We average across the
% orientation to get a population response from the neurons. The contrast %
% is used to examine the difference from the population firing rate and the 
% baseline background noise, indicated by an open parameter sigma. 


% (2) Experimental Design: Present stimulus of varying grating contrasts
% (from 0% to 75%) of 90 degree horizontal gradients. The individual neurons, at each
% contrast level, will fire more frequently when the preferred orientation
% is more visible (i.e. at higher contrast than at lower contrast).
% However, the firing rates of the neurons are normalized by an inhibitory
% layer of neuron which ensures that regardless of the input stimulus, the 
% population of neurons will fire and fire at a normalized or constrained
% rate that is insensitive to the change in input stimulus. 
%
% (3) Expected Results: We expect to see the effects of noralization using 
% this network architecture of having an
% inhibitory layer of neuron on the population reponses irrespective of the
% change in input visual stimulus. 
%
%%





%%

clear all
close all
clc

% Gaussian Tuning Curves: 
rmax= 52.14; % parameter from Dayan&Abbot 
sigma= 14.73; % parameter from Dayan&Abbot 
n=2; 
% Make Stimulus:
s = [-45:270]; % preferred orientation (degrees) 
c = [0;.12;.25;.50;.75]; 

% 7 Input Neurons: InpN 
% how do I create input neurons using cosine turning curves? 
PrefDir_Neuron = [0,22.5,45,90,112.5,135,180]; % input NEURON orientation angles in degrees

numInpN = length(PrefDir_Neuron);
numOutN = length(PrefDir_Neuron);

% allocate and initialize neurons
inpN = zeros(numInpN,length(s),length(c)); % input layer  
outN = zeros(numOutN,length(s),length(c)); % final output layer
outN1 = zeros(numOutN,length(s),length(c)); % output of first layer into inhibitory
outN2 = zeros(numOutN,length(s),length(c)); % output of inhibitory layer
inpINH = zeros(numInpN,length(s),length(c)); 
inpINH2 = zeros(numInpN,length(s),length(c));
inpOUT = zeros(numInpN,length(s),length(c));

% allocate and initialize 3 sets of weights: 
% 1. Excitatory weights from Input to Output neurons: 
    % the input to output weights are set to make horizaontal directions have
    % stronger weights projecting to the horizaontal output neuron and vice versa
    Vertical = [4]; 
%     w=ones(numInpN,1)*1.3;% initialized to .1
lower=0.5;
upper=1.5; 
w = lower + (upper-lower).*rand(numInpN,1); 
    w(Vertical) = 14; 
    wexc = repmat( w(:), [1,length(s),length(c)] );
    
% 2. InputNeurons to INHIBITORY neuron:
% winp = ones(numInpN,1)*-.05;
lower= -0.05;
upper= -0.1; 
winp = lower + (upper-lower).*rand(numInpN,1); 
wInpINH = repmat( winp(:),[1,length(s),length(c)] );

% 3. INHIBITORY neuron to OutputNeurons: 
% winh = ones(numInpN,1)*-1; % initialized to -1 
lower= -1;
upper= -1.1; 
winh = lower + (upper-lower).*rand(numInpN,1); 
wINHout = repmat( winh(:),[1,length(s),length(c)] ) ; 

% outN(1,:) = rand(1,length(s));    

for a=1:numInpN % neuron's PrefDir loop 
       for si=1:length(s) % loop through 361 degrees for inpN activity 
           for ci=1:length(c) % loop through the different contrast levels 
% calcuates the input derived from the stimulus 
inpN(a,si,ci) = rmax * exp( (-1/2)* ( ( s(si)-PrefDir_Neuron(a) ) ./ sigma ).^2 );
% excitatory input - output connection
inpOUT(a,si,ci) = wexc(a,si,ci) .* inpN(a,si,ci); 
% output activity from excitatory connection 
outN1(a,si,ci) = 1/(1+exp(( -1 .* inpOUT(a,si,ci))  ));
% inhibitory layer: 
inpINH(a,si,ci) = dot(wInpINH(a,si,ci), inpN(a,si,ci));
% inhibitory neuron output
outN2(a,si,ci) = 1/(1+exp(( -1 .* inpINH(a,si,ci))  ));
inpINH2(a,si,ci) = dot(wINHout(a,si,ci), outN2(a,si,ci));
% output from the inhibitory neuron 
outN22(a,si,ci) = 1/(1+exp(( -1 .* inpINH2(a,si,ci) )));
% take difference of excitatory and inhibitory outputs
outN(a,si,ci) = outN1(a,si,ci) - outN22(a,si,ci) ; 
% Normalize the output activities of the network by contrast % 
 norm = sigma^n + c(ci)^n;
outN(a,si,ci) = outN(a,si,ci)* c(ci) ./ norm; % ( c(ci) / (sigma + c(ci)) ); 
                
           end         
        end %   
end % end PrefDir loop  
    % mean firing rate of the population responses 
    mean_fr = sum(outN,2); 
        % interp1
    fr1=interp1(PrefDir_Neuron,mean_fr(:,1,1),[0:180],'cubic');
    figure;
    subplot(2,3,1);
    plot([0:180],fr1)
    title (['contrast=', num2str(c(1)),'%'])
    xlabel('degree')
    ylabel('activity')
%     ylim([0 2])
%     xlim([-90 270])
   
    subplot(2,3,2);
    fr2=interp1(PrefDir_Neuron,mean_fr(:,1,2),[0:180], 'cubic'); 
    plot([0:180],fr2)
        title (['contrast=', num2str(c(2)),'%'])
   
    xlabel('degree')
    ylabel('activity')
%     ylim([0 2])
%      xlim([-90 270])
     subplot(2,3,3);
      fr3=interp1(PrefDir_Neuron,mean_fr(:,1,3),[0:180], 'linear');
    plot([0:180],fr3)
        title (['contrast=', num2str(c(3)),'%'])
    xlabel('degree')
    ylabel('activity')
%     ylim([0 2])
%      xlim([-90 270])
     subplot(2,3,4);
    plot(PrefDir_Neuron,mean_fr(:,1,4),'*-')
        title (['contrast=', num2str(c(4)),'%'])
   
    xlabel('degree')
    ylabel('activity')
%     ylim([0 4]) 
%     xlim([-90 270])
     subplot(2,3,5);
    plot(PrefDir_Neuron,mean_fr(:,1,5),'*-')
        title (['contrast=', num2str(c(5)),'%']) 
    xlabel('degree')
    ylabel('activity')
% ylim([0 6])
%  xlim([-90 270])

no_norm = sum(outN1,2);
figure;
plot(PrefDir_Neuron,no_norm(:,1,1) )
    title (['WITHOUT NORMALIZATION '])
    xlabel('degree')
    ylabel('activity')
