%% Simulate stimuli (up direction 1 / down direction 0) 
% for the 100 trials at each of the 3 coherence levels
clc
clear all
close all

ntrials = 100; 
coherence = [0.8,3.2,12.8]% coherence level of the stimuli presented in the 100 trials 
lambdaup = [20;20;40]; % mean distribution of the firing rate responding to the UP preferred direction at each of the respective coherence level
lambdadown = [18;18;18];

for (n=1:length(coherence))
    
% stimuli = round(rand(ntrials,1)); % stimulus presented at each of the 100 trials
% 1 encode up direction, 0 encode down
upRate(n,:) = poissrnd(lambdaup(n),[ntrials,1])'; % poissand firing rate distribution for preferred PLUS direction 
downRate(n,:) = poissrnd(lambdadown(n), [ntrials,1])'; 

z = [0:max(upRate)]; % z threshold applied to discriminate the poissan firing rate distributions between UP and DOWN preferred directions   
% i.e. If r >= z, report "UP" otherwise "DOWN" direction preference 

pHit = zeros(length(z),length(coherence));
pFA =  zeros(length(z),length(coherence));

    for zi = 1:length(z)
        threshold = z(zi);
        pHit(zi,n) = sum(upRate(n,:)>threshold)/ntrials;
        pFA(zi,n) = sum(downRate(n,:)>threshold)/ntrials;
    end
    
figure;
clf
plot(pFA,pHit,'.-');

legend(num2str(coherence(n)),'Location','SouthEast');
set(gca,'XLim',[0,1]);
set(gca,'YLim',[0,1]);

title('ROC Curve for z threshold at each motion coherence level','FontSize',12,'FontWeight','bold')
xlabel('\alpha (false alarm rate)','FontSize',12,'FontWeight','bold') % x-axis label
ylabel('\beta (hit rate)','FontSize',12,'FontWeight','bold') % y-axis label  
set(gca,'XTick',0:.2:1);
set(gca,'YTick',0:.2:1);
hold on
plot([0,1],[0,1],'k-');
axis square
   
end










%%




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








