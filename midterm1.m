%% midterm #1 Divisive Normalization   author: Colleen Chen
% (1) Implementation: 
%
% (2) Experimental Design: Present stimulus of grating contrasts of varying
% lengh
%
% (3) Expected Results: 
%
%


clear all
close all
clc

% Gaussian Tuning Curves: 
rmax= 52.14; % parameter from Dayan&Abbot 
sigma= 14.73; % parameter from Dayan&Abbot 
% Make Stimulus:
s = [-90:270]; % preferred orientation (degrees) 
c = [0;.12;.25;.50;.75]; 

% 7 Input Neurons: InpN 
% how do I create input neurons using cosine turning curves? 
PrefDir_Neuron = [0,22.5,45,90,112.5,135,180]; % input NEURON orientation angles in degrees

numInpN = length(PrefDir_Neuron);
numOutN = length(PrefDir_Neuron);

% allocate and initialize neurons
inpN = zeros(numInpN,length(s),length(c)); %  
outN = zeros(numOutN,length(s),length(c));
inpINH = zeros(numInpN,length(s),length(c)); 


% allocate and initialize 3 sets of weights: 
% 1. Excitatory weights from Input to Output neurons: 
    % the input to output weights are set to make horizaontal directions have
    % stronger weights projecting to the horizaontal output neuron and vice versa
    Vertical = 4; 
    w=ones(numInpN,1)*.1;% initialized to .1
    w(Vertical) = 3; 
    wexc = repmat( w(:), [1,length(s),length(c)] );
    

% 2. InputNeurons to INHIBITORY neuron:
winp = ones(numInpN,1)*.5;
wInpINH = repmat( winp(:),[1,length(s),length(c)] );

% 3. INHIBITORY neuron to OutputNeurons: 
winh=ones(numInpN,1)*-1; % initialized to -1 
wINHout = repmat( winh(:),[1,length(s),length(c)] ) ; 

% outN(1,:) = rand(1,length(s));    

for a=1:numInpN % neuron's PrefDir loop 
       for si=1:length(s) % loop through 361 degrees for inpN activity 
           for ci=1:length(c) % loop through the different contrast levels 

                inpN(a,si,ci) = rmax * c(ci) * exp( (-1/2)* ( ( s(si)-PrefDir_Neuron(a) ) ./ sigma ).^2 );
                inpINH(a,si,ci) = wexc(a,si,ci) .* inpN(a,si,ci); 
                outN(a,si,ci) = inpINH(a,si,ci) .* ( c(ci)/(sigma+c(ci)) );
                
           end         
        end % 
       
end % end PrefDir loop  
   
mean_fr = sum(outN,2); 
    figure;
    plot(PrefDir_Neuron,mean_fr(:,1,1))
    hold on 
    plot(PrefDir_Neuron,mean_fr(:,1,2))
    plot(PrefDir_Neuron,mean_fr(:,1,3))
    plot(PrefDir_Neuron,mean_fr(:,1,4))
    plot(PrefDir_Neuron,mean_fr(:,1,5))
    
    
    subplot(2,3,1);
    plot(s,outN(1,:,1))
    hold on 
    plot(s,outN(2,:,1))
    plot(s,outN(3,:,1))
    plot(s,outN(4,:,1))
    plot(s,outN(5,:,1))
    plot(s,outN(6,:,1))
    plot(s,outN(7,:,1))
    
    subplot(2,3,2);
    plot(s,outN(1,:,2))
    hold on 
    plot(s,outN(2,:,2))
    plot(s,outN(3,:,2))
    plot(s,outN(4,:,2))
    plot(s,outN(5,:,2))
    plot(s,outN(6,:,2))
    plot(s,outN(7,:,2))
    
    subplot(2,3,3);
    plot(s,outN(1,:,3))
    hold on
    plot(s,outN(2,:,3))
    plot(s,outN(3,:,3))
    plot(s,outN(4,:,3))
    plot(s,outN(5,:,3))
    plot(s,outN(6,:,3))
    plot(s,outN(7,:,3))
    
    subplot(2,3,4);
    plot(s,outN(1,:,4))
    hold on 
    plot(s,outN(2,:,4))
    plot(s,outN(3,:,4))
    plot(s,outN(4,:,4))
    plot(s,outN(5,:,4))
    plot(s,outN(6,:,4))
    plot(s,outN(7,:,4))


    subplot(2,3,5);
    plot(s,outN(1,:,5))
    hold on 
    plot(s,outN(2,:,5))
    plot(s,outN(3,:,5))
    plot(s,outN(4,:,5))
    plot(s,outN(5,:,5))
    plot(s,outN(6,:,5))
    plot(s,outN(7,:,5))
    
    
%%
% clear all
% close all
% clc
% 
% % Gaussian Tuning Curves: 
% rmax= 52.14; % parameter from Dayan&Abbot 
% sigma= 14.73; % parameter from Dayan&Abbot 
% % Make Stimulus:
% s = [-90:270]; % preferred orientation (degrees) 
% c = [0;.12;.25;.50;.75]; 
% 
% % 7 Input Neurons: InpN 
% % how do I create input neurons using cosine turning curves? 
% PrefDir_Neuron = [0,22.5,45,90,112.5,135,180]; % input NEURON orientation angles in degrees
% 
% numInpN = length(PrefDir_Neuron);
% numOutN = length(PrefDir_Neuron);
% 
% % allocate and initialize neurons
% inpN = zeros(numInpN,length(s),length(c)); %  
% outN = zeros(numOutN,length(s),length(c));
% inpINH = zeros(numInpN,length(s),length(c)); 
% 
% 
% % allocate and initialize 3 sets of weights: 
% % 1. Excitatory weights from Input to Output neurons: 
%     % the input to output weights are set to make horizaontal directions have
%     % stronger weights projecting to the horizaontal output neuron and vice versa
%     Vertical = 4; 
%     w=ones(numInpN,1)*.1;% initialized to .1
%     w(Vertical) = 3; 
%     wexc = repmat( w(:), [1,length(s),length(c)] );
%     
% 
% % 2. InputNeurons to INHIBITORY neuron:
% winp = ones(numInpN,1)*.5;
% wInpINH = repmat( winp(:),[1,length(s),length(c)] );
% 
% % 3. INHIBITORY neuron to OutputNeurons: 
% winh=ones(numInpN,1)*-1; % initialized to -1 
% wINHout = repmat( winh(:),[1,length(s),length(c)] ) ; 
% 
% % outN(1,:) = rand(1,length(s));    
% 
% for a=1:numInpN % neuron's PrefDir loop 
%        for si=1:length(s) % loop through 361 degrees for inpN activity 
%            for ci=1:length(c) % loop through the different contrast levels 
% 
%                 inpN(a,si,ci) = rmax * c(ci) * exp( (-1/2)* ( ( s(si)-PrefDir_Neuron(a) ) ./ sigma ).^2 );
%                 inpINH(a,si,ci) = wexc(a,si,ci) .* inpN(a,si,ci); 
%                 outN(a,si,ci) = inpINH(a,si,ci) .* ( c(ci)/(sigma+c(ci)) );
%                 
%            end         
%         end % 
%        
% end % end PrefDir loop  
%    
% mean_fr = sum(outN,2); 
%     figure;
%     plot(PrefDir_Neuron,mean_fr(:,1,1))
%     hold on 
%     plot(PrefDir_Neuron,mean_fr(:,1,2))
%     plot(PrefDir_Neuron,mean_fr(:,1,3))
%     plot(PrefDir_Neuron,mean_fr(:,1,4))
%     plot(PrefDir_Neuron,mean_fr(:,1,5))
%     
%     
%     subplot(2,3,1);
%     plot(s,outN(1,:,1))
%     hold on 
%     plot(s,outN(2,:,1))
%     plot(s,outN(3,:,1))
%     plot(s,outN(4,:,1))
%     plot(s,outN(5,:,1))
%     plot(s,outN(6,:,1))
%     plot(s,outN(7,:,1))
%     
%     subplot(2,3,2);
%     plot(s,outN(1,:,2))
%     hold on 
%     plot(s,outN(2,:,2))
%     plot(s,outN(3,:,2))
%     plot(s,outN(4,:,2))
%     plot(s,outN(5,:,2))
%     plot(s,outN(6,:,2))
%     plot(s,outN(7,:,2))
%     
%     subplot(2,3,3);
%     plot(s,outN(1,:,3))
%     hold on
%     plot(s,outN(2,:,3))
%     plot(s,outN(3,:,3))
%     plot(s,outN(4,:,3))
%     plot(s,outN(5,:,3))
%     plot(s,outN(6,:,3))
%     plot(s,outN(7,:,3))
%     
%     subplot(2,3,4);
%     plot(s,outN(1,:,4))
%     hold on 
%     plot(s,outN(2,:,4))
%     plot(s,outN(3,:,4))
%     plot(s,outN(4,:,4))
%     plot(s,outN(5,:,4))
%     plot(s,outN(6,:,4))
%     plot(s,outN(7,:,4))
% 
% 
%     subplot(2,3,5);
%     plot(s,outN(1,:,5))
%     hold on 
%     plot(s,outN(2,:,5))
%     plot(s,outN(3,:,5))
%     plot(s,outN(4,:,5))
%     plot(s,outN(5,:,5))
%     plot(s,outN(6,:,5))
%     plot(s,outN(7,:,5))
%     
%     