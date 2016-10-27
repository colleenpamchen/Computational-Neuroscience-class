%% Homework3.m 
time = [0:100]; 
rnd = a + (b-a).*rand( length() ,1); 
sigmoid_activation_function = 1 / ( 1 + exp(-1 * gain * input) )



for t=0:100 % loop through t 
 
% Excitatory neurons E1 & E2 
E1(t+1) = N_extrinsic1 * W_extrinsic + E2(t) * W_intrinsic - E2(t) * W_inhibitory + rnd ; 
E2(t+1) = N_extrinsic2 * W_extrinsic + E1(t) * W_intrinsic - E1(t) * W_inhibitory + rnd ; 

end


plot(time,E1,'k:')
hold on 
plot(time,E2,'r-')

