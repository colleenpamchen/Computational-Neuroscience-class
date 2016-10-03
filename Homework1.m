
%% Psy 268A Homework 1
% Vary the width of the Gaussian Tuning Curve 
clc
clear all
close all

s = -40:1:40;
rmax= 52.14;
smax= 0;
sig= [14.73, (14^2), (1/14)];

for i=1:length(sig) 
y = rmax * exp( (-1/2)* ( (s - smax) / sig(i) ).^2 );
figure
plot(s,y)
title('Gaussian tuning curve of avg. firing rate','FontSize',12,'FontWeight','bold')
xlabel('s: orientation angle in degrees','FontSize',12,'FontWeight','bold') 
ylabel('f(Hz)','FontSize',12,'FontWeight','bold') 
legend( num2str(sig(i)) ,'Location','SouthEast');
end

%% Vary the gain rate of Sigmoid Tuning Curve 
clc
clear all
close all

s = -1:.01:1; % retinal disparity 
rmax = 36.03; 
shalf = 0.036; % disparity that produces a firing rate half as big as the maximum value rmax
deltas = [0.029, 0.06, 0.12, 0.2]; % rate of incrase 
y = zeros(1,length(s));

% signmoid (logit) function
for (n=1:length(deltas))

for (i=1:length(s))
y(i) = rmax / (1+exp( ( shalf - s(i) ) / deltas(n) ) );
end

figure
plot(s,y)
title('Sigmoid turning curve of retinal disparity','FontSize',12,'FontWeight','bold')
xlabel('s: retinal disparity in degrees','FontSize',12,'FontWeight','bold') 
ylabel('f(Hz)','FontSize',12,'FontWeight','bold') 
legend( num2str(deltas(n)) ,'Location','SouthEast');

end


%% ROC curve for simulated discrimination analysis 
% simulated 100 trials at each of the 3 coherence levels
clc
clear all
close all

ntrials = 100; 
coherence = [0.8,3.2,12.8] % coherence level of the stimuli presented in the 100 trials 
lambdaup = [20;20;40]; % mean distribution of the firing rate responding to the preferred direction at each of the respective coherence level
lambdadown = [18;18;18];

for (n=1:length(coherence))
    upRate(n,:) = poissrnd(lambdaup(n),[ntrials,1])'; % poissand firing rate distribution for preferred direction (UP) 
    downRate(n,:) = poissrnd(lambdadown(n), [ntrials,1])'; % firing rate distribution for Nonpreferred direction (DOWN)
    z = [0:max(upRate)]; % z threshold applied to discriminate the poissan firing rate distributions between preferred and Nonpreferred directions   
    % i.e. If r >= z, report "UP" otherwise "DOWN" directional preference 

    pHit = zeros(length(z),length(coherence)); % hit rate 
    pFA =  zeros(length(z),length(coherence)); % false alarm rate 

        for zi = 1:length(z)
            threshold = z(zi);
            pHit(zi,n) = sum(upRate(n,:)>threshold)/ntrials; % probability of hit rate
            pFA(zi,n) = sum(downRate(n,:)>threshold)/ntrials; % probability of false alarm rate
        end

    figure;
    clf
    plot(pFA,pHit,'.-');

    legend(num2str(coherence(n)),'Location','SouthEast');
    set(gca,'XLim',[0,1]);
    set(gca,'YLim',[0,1]);

    title('ROC Curve for z threshold at each motion coherence level','FontSize',12,'FontWeight','bold')
    xlabel('\alpha (false alarm rate)','FontSize',12,'FontWeight','bold') 
    ylabel('\beta (hit rate)','FontSize',12,'FontWeight','bold')   
    set(gca,'XTick',0:.2:1);
    set(gca,'YTick',0:.2:1);
    hold on
    plot([0,1],[0,1],'k-');
    axis square
   
end












