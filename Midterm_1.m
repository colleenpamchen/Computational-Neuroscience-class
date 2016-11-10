%% Midterm #1 Divisive Normalization
% (1) Implementation: 
%
% (2) Experimental Design: 
%
% (3) Expected Results: 
%
%

function cc_midterm1

sigma= 14.73;

% Make Stimulus:
s = [-90:270]; % preferred orientation (degrees) 

contrast = [0;12;25;50;75]; % ??????????? HOW to use contrast with pd ?? 
Time=100; 

% 7 Input Neurons: InpN 
% how do I create input neurons using cosine turning curves? 
Npd = [0,22.5,45,90,112.5,135,180]; % input NEURON orientation angles in degrees

numInp = length(Npd);
numOut = length(Npd);

% allocate and initialize neurons
inpN = zeros(numInp,length(s)); % inpN =zeros(length(pd),length(angle)); 
outN = zeros(Time,length(s));

% allocate and initialize 3 sets of weights: 
% 1. Excitatory weights from Input to Output neurons: 

wext=[1,1,1,1,1,1,1];% initialized to 1
wInpOutExc = repmat(wext(:),1,length(s));

    % the input to output weights are set to make horizaontal directions have
    % stronger weights projecting to the horizaontal output neuron and vice versa
    Vertical = 4; 
    wInpOutExc(Vertical,:) = 10; 

% 2. InputNeurons to INHIBITORY neuron:
winp = [.5,.5,.5,.5,.5,.5,.5];
wInpINH = repmat(winp(:),1,length(s));

% 3. INHIBITORY neuron to OutputNeurons: 
winh=[-1,-1,-1,-1,-1,-1,-1]; % initialized to -1 
wINHout = repmat(winh(:),1,length(s)) ; 

outN(1,:) = rand(1,length(s));    

for a=1:numInp % neuron's ANGLE loop 
       for p=1:length(s) % pd loop for inpN activity 
            inpN(a,p) = cosTune(s(p), Npd(a));
            
       end
       
% figure;
% plot(Npd,inpN(a,:))
% hold on
% contrast(a,:) = inpN(a,:) .* (0.5 / (sigma + 0.5) );
% plot(Npd,contrast(a,:)); 
% hold off 

       for t=2:Time
    % set the rates of the input neurons
     
        % cycle through the output neurons
        for i=1:numOut % loop through OutputN indexed 'i' 
            % initialize SynIn   
              SynIn = 0 ; 
            for j = 1:numInp % loop through InputN indexed 'j' 
                % (1) 1-1 Input -> Output EXC : SynIn = SynIn + InpN .* wInpOutExc
%                 SynIn =  SynIn + inpN(j) .* wInpOutExc(j,i) ;
                SynIn = inpN(j,:) .* wInpOutExc(j,:) ;
                % (2) all-1 Input -> INH : SynIn = SynIn + InpN .* wInpINH
%                 SynIn = SynIn + inpN(j) .* wInpINH(j) ;  
                SynIn = SynIn + inpN(j,:) .*  wInpINH(j,:) ;
                 % (3) 1-all INH -> Output : SynIn = SynIn + InpN .* wINHout
                SynIn = SynIn + dot(outN(t-1,j), wINHout(j)) ;   
                    
            end

%             outN(t,i) = sigmoidN(SynIn);
            outN(t,:) = SynIn;
             
            outfr(i) = sum(outN,i); 
        end % outN 
        
    end % end of time loop
%     outfr(a) = sum(outN(:,a));
%     figure
%     plot(Npd,outN(:,i),'-*')
end % end neuron ANGLE loop 
    
    figure
    plot(Npd,outfr,'-*')
    

end % end of user-defined function

function r = cosTune(pd, s)
% r = abs(cos(s-pd));
% Gaussian Tuning Curves: 
rmax= 52.14; % parameter from Dayan&Abbot 
sig= 14.73; % parameter from Dayan&Abbot 
r = rmax * exp( (-1/2)* ( (pd-s)/sig ).^2 );
end
% Sigmoid Activation
%   Input:  I - synaptic input
%           g - gain
%   Output: firing rate of sigmoid function
% function r = sigmoidN (I) 
% g=1; 
% r = 36.03/(1+exp((-I)/.029  ));
% end

