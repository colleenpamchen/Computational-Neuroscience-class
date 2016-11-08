function MacHW3

itype = 2;

%%% first layer %%%
% cosine tuning parameters
degs = 1:.1:360; tw=1; r0 = 0; rmax = 1.5;
% preferred orientations 
prefs = [45,90,135,180];
% initialize neuron response matrix;
n = nan(length(degs),length(prefs));
% generate firing rates;
for pref = prefs
 n(:,(prefs==pref)) = costu(r0,rmax,pref,tw,degs);
end

%%% second layer %%%
wv = [0;1;0;0]; % weights for vertical neuron
wh = [0;0;0;1]; % weights for horizontal neuron
rvn = sig(n*wv);% rate for veritcal neuron, no inhibition
rhn = sig(n*wh);% rate for horizontal neuron, no inhibition

%%% different ways to express inhibition...
switch itype
 case 1
  ivn= max(rvn./(rhn./max(rhn)),0); % vert response after inhibition
  ihn= max(rhn./(rvn./max(rvn)),0); % horz response after inhibition
 case 2 
  ivn = max(rvn-rhn,0);
  ihn = max(rhn-rvn,0);
 case 3
  ivn = sig(n*wv-rhn);
  ihn = sig(n*wh-rvn);
end

%%% plotting %%%
figure;
% Raw neurons;
subplot(3,1,1); hold on;
plot(degs,n)
h = legend('45 deg','90 deg','135 deg','180 deg');
set(h,'location','best')
title('neuronal responses')
xlabel('degrees')
ylabel('rate (arbitrary units)')
% raw second layer
subplot(3,1,2); hold on;
maxr = max([rvn(:);rhn(:)]);
plot([90,270],[maxr,maxr],'k*');
plot(degs,rvn)
plot([0,180],[maxr,maxr],'ko');
plot(degs,rhn)
h=legend('vertical orientation','vartical response',...
         'horizontal orientation','horizontal response');
set(h,'location','best')
xlabel('degrees')
ylabel('response (arbitrarty units)')
title('Vertical and horizontal neurons w/o inhibition')
% raw second layer;
subplot(3,1,3); hold on;
maxi = max([ivn(:);ihn(:)]);
plot([90,270],[maxi,maxi],'k*'); 
plot(degs,ivn);
plot([0,180],[maxi,maxi],'ko');
plot(degs,ihn);
h = legend('vertical orientation','vertical neuron',...
           'horizontal orientation','horizontal neuron');
set(h,'location','best');
ylabel('firing rate (arbitrary units)')
xlabel('angle (degrees)')
title('Vertical and horizontal responses with inhibition')

end

% cosine tuning curves, takes degree args
function rates = costu(r0,rmax,pref,tw,degs)
rates = max(0, r0+(rmax-r0)*cosd(degs-pref).^tw);
end

% sigmoid function
function s = sig(x)
s = 1./(1+exp(-x));
end



