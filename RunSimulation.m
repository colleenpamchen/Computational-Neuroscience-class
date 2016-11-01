function RunSimulation
% Uniformly sample the simulation parameter space.  
gain = 1:10;
weights1 = [.2,.4,.8,1.6,2.4,3.2,4.8,6.4,9.6,12.8];
weights2 = [.2,.4,.8,1.6,2.4,3.2,4.8,6.4,9.6,12.8];
activityDiff = [.01,.05,.075,.1,.125,.15,.175,.2,.225,.25];
[C1,C2] = meshgrid(activityDiff,gain);
gainConds = [C1(:),C2(:)];
[C1,C2] = meshgrid(activityDiff,weights1);
w1Conds = [C1(:),C2(:)];
[C1,C2] = meshgrid(activityDiff,weights2);
w2Conds = [C1(:),C2(:)];

nSims = length(gainConds(:,1));
totalDur = 100; % milliseconds
nReps = 100;
ds = .2;  %
Amax = ones(2,1);
halfMax = zeros(2,1);

halfTime = NaN(nReps,nSims,3);
endTime = NaN(nReps,nSims,3);
excitation = NaN(2,totalDur,nSims,3);
firingRate = NaN(2,totalDur,nSims,3);
successes = NaN(nSims,2,3);
for sim = 1:nSims
    afferentSignal = .5*ones(2,1);
    afferentSignal(2) =  afferentSignal(2) + w1Conds(sim,1);
    gainVec = ones(2,1);
    W(1) = w1Conds(sim,2);
    W(2) = .1;
    W(3) = w1Conds(sim,2);
    [successes(sim,:,1),halfTime(:,sim,1),endTime(:,sim,1),excitation(:,:,sim,1),firingRate(:,:,sim,1)] = runTrial(nReps,totalDur,gainVec,Amax,halfMax,W,ds,afferentSignal);

    afferentSignal = .5*ones(2,1);
    afferentSignal(2) =  afferentSignal(2) + w2Conds(sim,1);
    gainVec = ones(2,1);
    W(1) = .1;
    W(2) = w2Conds(sim,2);
    W(3) = .1;
    [successes(sim,:,2),halfTime(:,sim,2),endTime(:,sim,2),excitation(:,:,sim,2),firingRate(:,:,sim,2)] = runTrial(nReps,totalDur,gainVec,Amax,halfMax,W,ds,afferentSignal);

    afferentSignal = .5*ones(2,1);
    afferentSignal(2) =  afferentSignal(2) + gainConds(sim,1);
    gainVec = gainConds(sim,2)*ones(2,1);
    W(1) = .1;
    W(2) = .1;
    W(3) = .1;
    [successes(sim,:,3),halfTime(:,sim,3),endTime(:,sim,3),excitation(:,:,sim,3),firingRate(:,:,sim,3)] = runTrial(nReps,totalDur,gainVec,Amax,halfMax,W,ds,afferentSignal);
end
% add ytick labels.
fig1A = zeros(length(activityDiff),length(weights1)); % Figure 4a, left
fig1B = zeros(length(activityDiff),length(weights1)); %figure 4a, right
for w = 1:length(weights1)
   dltW = weights1(w); 
   for a = 1:length(activityDiff)
      dltA = activityDiff(a);
      indx = w1Conds(:,1) == dltA & w1Conds(:,2) == dltW;
      fig1A(a,w) = successes(indx,1,1);
      fig1B(a,w) = successes(indx,2,1);
   end
end
figure;
subplot(1,2,1);
imagesc(fig1A);
colorbar
title('Modulate Ext/Inh Weights: Time Step 50')
xlabel('Extrinsic and Inhibitory Weights')
ylabel('Activity Difference')
set(gca,'YTickLabel',activityDiff);
subplot(1,2,2);
imagesc(fig1B);
colorbar
title('Modulate Ext/Inh Weights: Time Step 100')
xlabel('Extrinsic and Inhibitory Weights')
ylabel('Activity Difference')
set(gca,'YTickLabel',activityDiff);

fig2A = zeros(length(activityDiff),length(weights2)); % Figure 4b, left
fig2B = zeros(length(activityDiff),length(weights2)); %figure 4b, right
for w = 1:length(weights1)
   dltW = weights2(w); 
   for a = 1:length(activityDiff)
      dltA = activityDiff(a);
      indx = w2Conds(:,1) == dltA & w2Conds(:,2) == dltW;
      fig2A(a,w) = successes(indx,1,2);
      fig2B(a,w) = successes(indx,2,2);
   end
end
figure;
subplot(1,2,1);
imagesc(fig2A);
colorbar
title('Modulate Intrinsic: Time Step 50')
xlabel('Intrinsic Weights')
ylabel('Activity Difference')
set(gca,'YTickLabel',activityDiff);
subplot(1,2,2);
imagesc(fig2B);
colorbar
title('Modulate Intrinsic: Time Step 100')
xlabel('Intrinsic Weights')
ylabel('Activity Difference')
set(gca,'YTickLabel',activityDiff);

fig3A = zeros(length(activityDiff),length(gain)); % Figure 4c, left
fig3B = zeros(length(activityDiff),length(gain)); %figure 4c, right
for w = 1:length(weights1)
   dltW = gain(w); 
   for a = 1:length(activityDiff)
      dltA = activityDiff(a);
      indx = gainConds(:,1) == dltA & gainConds(:,2) == dltW;
      fig3A(a,w) = successes(indx,1,2);
      fig3B(a,w) = successes(indx,2,2);
   end
end
figure;
subplot(1,2,1);
imagesc(fig3A);
colorbar
title('Modulate Gain: Time Step 50')
xlabel('Gain')
ylabel('Activity Difference')
set(gca,'YTickLabel',activityDiff);
subplot(1,2,2);
imagesc(fig3B);
colorbar
title('Modulate Gain: Time Step 100')
xlabel('Gain')
ylabel('Activity Difference')
set(gca,'YTickLabel',activityDiff);

end

%% subfunction
function [success,halfTime,endTime,outdV,outRate] = runTrial(nReps,totalDur,gainVec,Amax,halfMax,W,ds,afferentSignal)
% within function fixed parameters
nExtUnits = 2;
nInhUnits = 2; 
%preallocate outputs for speed.  
firingRate = zeros(2,totalDur,nReps);
halfTime = NaN(nReps,1);
endTime = NaN(nReps,1);
dV = NaN(2,totalDur,nReps);
% internal function: Activation 
activationFunction = @(input,gain,Amax,halfMax)...
Amax./(1+exp(-gain.*(halfMax-input)));
% repeat trial nRep times.  
for rep = 1:nReps
        E = zeros(2,totalDur);
        rate = zeros(2,totalDur);
        Wext = .1*eye(nExtUnits);
        Wint = .1*flipud(eye(nExtUnits));
        Winh = .1*flipud(eye(nInhUnits));
        gain = ones(2,1);
for t = 1:totalDur
    if t == totalDur/2 
         % phasic mode
        Wext = W(1)*eye(nExtUnits);
        Wint = W(2)*flipud(eye(nExtUnits));
        Winh = W(3)*flipud(eye(nInhUnits));
        gain = gainVec;
    end
    input = ds*afferentSignal;  
    s = 0;
    while s < 1
        if s == 0
           tmpE = Wext*input + Wint*E(:,t) - Winh*E(:,t) ;
        else
           tmpE = Wext*input + Wint*tmpE - Winh*tmpE;
        end
       s = s + ds;
    end
    E(:,t+1) = tmpE + 2*rand(2,1)-1;
    rate(:,t+1) = rate(:,t) + activationFunction(E(:,t+1),gain,Amax,halfMax);
end
dV(:,:,rep) = E(:,1:totalDur);
firingRate(:,:,rep) = rate(:,1:totalDur);
halfTime(rep) = rate(2,totalDur/2) >= rate(1,totalDur/2);
endTime(rep) = rate(2,end) >= rate(1,end);
end
outdV = mean(dV,3);
outRate = mean(firingRate,3);
success(1) = sum(halfTime);
success(2) = sum(endTime);
end