%% Simulate stimuli (up direction 1 / down direction 0) 
% for the 100 trials at each of the 3 coherence levels

ntrials = 100; 

coherence = [0.8; 3.2; 12.8]; % coherence level of the stimuli presented in the 100 trials 
lambda = [20; 20; 40]; % mean distribution of the firing rate responding to the UP preferred direction at each of the respective coherence level

z = [1:100]'; % z threshold applied to discriminate the poissan firing rate distributions between UP and DOWN preferred directions   
% i.e. If r >= z, report "UP" otherwise "DOWN" direction preference 

for (n=1:length(coherence))
    
stimuli = round(rand(ntrials,1)); % stimulus presented at each of the 100 trials 
r = poissrnd(lambda(n),[ntrials,1]); % poissand firing rate distribution for preferred PLUS direction 

for (i=1:length(stimuli))
predicted(i) = ( r(i) >= z(i) ); % Predicted behavioral response based on applying r, firing rate data, to a z threshold. 
% UP direction is encoded 1; DOWN is 0 
hit(i) = logical(predicted(i) == stimuli(i)); % hit 
% fa = logical( find(stimuli == 0 & predicted ==1)); % false alarm
fa(i) = logical( stimuli(i) == 0 & predicted ==1 );

% end


beta(i) = sum(hit)/ntrials; 
alpha(i) = sum(fa)/ntrials; 


figure(n); 
title('ROC Curve for z threshold at each motion coherence level','FontSize',12,'FontWeight','bold')
axis([0,1,0,1])
xlabel('\alpha (false alarm rate)','FontSize',12,'FontWeight','bold') % x-axis label
ylabel('\beta (hit rate)','FontSize',12,'FontWeight','bold') % y-axis label
legend('% coherence level')

hold on ;
plot(alpha(i), beta(i),'x');


end



end


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








