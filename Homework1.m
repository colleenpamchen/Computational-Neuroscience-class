%% Simulate stimuli (up direction 1 / down direction 0) 
% for the 100 trials at each of the 3 coherence levels

ntrials = 100; 

coherence = [0.8; 3.2; 12.8];
lambda1 = 40; % for coherence 12.8%

% stimuli = struct (coherence, round(rand(1,100)) );
stimuli = round(rand(ntrials,1));  

r=poissrnd(lambda1,[ntrials,1]); 

z =25; 

predicted = ( r >=z ); 
hit = (predicted == stimuli); % hit 
fa = logical( find(stimuli == 0 & predicted ==1)); % false alarm

beta = sum(hit)/ntrials; 
alpha = length(fa)/ntrials; 


figure
plot(alpha, beta,'x');

title('ROC Curve for z threshold at each motion coherence level','FontSize',12,'FontWeight','bold')
axis([0,1,0,1])
xlabel('\alpha (false alarm rate)','FontSize',12,'FontWeight','bold') % x-axis label
ylabel('\beta (hit rate)','FontSize',12,'FontWeight','bold') % y-axis label
legend('% coherence level','coherence level','FontSize',12,'FontWeight','bold')
 


%% for each of the 100 stimuli, at each coherence level, apply the z threshold to the 
% averaged firing rate r (which comes from the neural data, ranged 0-60) to predict behavioral response 
% for directional prefernce (+/-). 




%% calculate the alpha (hit rate) and beta (false alarm) rates as datapoints over 100 trials




%% Psy 268A Homework 1

s = -40:1:40;
rmax= 52.14;
smax= 0;
sig= 14.73;

y = rmax * exp( (-1/2)* (s - smax / sig).^2 );

figure
plot(s,y)
title('Gaussian turning curve of avg. firing rate')
xlabel('s: orientation angle in degrees') % x-axis label
ylabel('f(Hz)') % y-axis label

