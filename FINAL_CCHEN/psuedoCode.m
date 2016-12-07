%% STDP and homeostasis Algorithm:
% Initialize variables:
Izhikevich 
STDP
Homeostasis 

% Generate Poisson spike times for each of the 100 input neurons

loop: sec= 1:1000 % loop through 1000 seconds of simulation

	loop: t=1:1000 % loop through ms 
    
		If ( Poisson spiked) 
            keep index of input(i) that fired ; 
			Generate next spike time for those input(i) that fired ;
			Update: I = s(fired) ;
			Update voltage: v=v+(0.04*v^2+5*v+140-u+ sum(I(fired)) ; 
			If (v >=30) % check to see if output neuron filed 
				R+=1; % keep count and increment firing rate
				v=c; % reset v 
				% apply STDP (LTP+LTD) with Homeostasis term K 
				K = R / (T.*(1+abs(1-R/Rtarget)*gamma));
           	 	dsdt = dsdt+ (alpha .* S1 .* (1 - R/Rtarget) + 1*(LTP1(fired)))*K;
			End
            dsdt = dsdt+ (alpha .* S1 .* (1 - R/Rtarget) + 1*(LTD1(fired)))*K;
            % exponentially decay LTP & LTD
	End
End
% update weight S
S = S + dsdt; 
S = max(wmin, S);

Plot firing rate as a function of time
Plot weights S of each 100 synapse 
	




		


