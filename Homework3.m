%% Homework3.m 
close all 
clear all

time = [1:101]; 
c=-1;
b=1; 
rnd = c + (b-c).*rand( length(time) ,1); 
gain = 1; % TUNE THIS parameter in the next step  
sig = @(input) 1./( 1 + exp(-1 * gain * input ));  

E1 = zeros(length(time) ,1);
E2 = zeros(length(time) ,1);
rate1 = zeros(length(time),1);
rate2 = zeros(length(time),1);

W_extrinsic = 0.1; 
W_intrinsic = 0.1; 
W_inhibitory = 0.1; 
% first iteration, when 
N_extrinsic1 = 0.5;
N_extrinsic2 = 0.6;
% initial condition
E1(1)=0.1;
E2(1)=0.9;
t=1;
E1(t+1) = N_extrinsic1 * W_extrinsic + E2(t) * W_intrinsic - E2(t) * W_inhibitory + rnd(t) ; 
rate1(t+1) =  sig(E1(t+1));

E2(t+1) = N_extrinsic2 * W_extrinsic + E1(t) * W_intrinsic - E1(t) * W_inhibitory + rnd(t) ; 
rate2(t+1) = sig(E2(t+1));

for t=2:floor(length(time)/2) % FIRST 50 simulation cycles: 
E1(t+1) = N_extrinsic1 * W_extrinsic + E2(t) * W_intrinsic - E2(t) * W_inhibitory + rnd(t) ; 
rate1(t+1) = rate1(t) + sig(E1(t+1));

E2(t+1) = N_extrinsic2 * W_extrinsic + E1(t) * W_intrinsic - E1(t) * W_inhibitory + rnd(t) ; 
rate2(t+1) = rate2(t) + sig(E2(t+1));

end

% LAST 50 simulation cycles:
W_extrinsic = 2.5;
W_inhibitory = 2.5; 

for t=ceil(length(time)/2):length(time)   
E1(t+1) = N_extrinsic1 * W_extrinsic + E2(t) * W_intrinsic - E2(t) * W_inhibitory + rnd(t) ; 
rate1(t+1) = rate1(t) + sig(E1(t+1));

E2(t+1) = N_extrinsic2 * W_extrinsic + E1(t) * W_intrinsic - E1(t) * W_inhibitory + rnd(t) ;
rate2(t+1) = rate2(t) + sig(E2(t+1));

end
 figure
        plot(time,rate1(2:102),'k')
        xlim([0 100])
        title('Activities of neurons E1 and E2' )
        xlabel('time')
        ylabel('firing rate')
        hold on 
        plot(time,rate2(2:102),'r')
        legend('E1','E2')
        hold off
        
% keep track of firing rates of E2 > E1 at HALFWAY and END of each trial
% store the ratio of E2 > E1 from the 100 runs 
% "activity difference" refers to: EXT2 - EXT1 

a = [0.05:0.05:0.25]; % activity difference, calculated from EXT2-EXT1 
W = [0.4, 1.6, 2.5, 3.6, 6.4, 10];
W_fixed = 0.1;
M50 = zeros(length(a),length(W));
M100 = zeros(length(a),length(W));

t=1; % reset t

for i = 1:length(a)
    
    for j = 1:length(W)
            tmp1= zeros(1,100);
            tmp2= zeros(1,100);
        for trial = 1:100

            
            for t = 2:length(time)
                N_extrinsic1=0.5; 
                N_extrinsic2 = a(i) + N_extrinsic1; 
                
                E1(t+1) = N_extrinsic1 * W_fixed + E2(t) * W_fixed - E2(t) * W_fixed + (c + (b-c).*rand(1))  ; % E1
                rate1(t+1) = rate1(t) + sig(E1(t+1));

                E2(t+1) = N_extrinsic2 * W_fixed + E1(t) * W_fixed - E1(t) * W_fixed + (c + (b-c).*rand(1)) ; % E2
                rate2(t+1) = rate2(t) + sig(E2(t+1));
                             
                if t > 52  % midtrial, 50th time step     
                E1(t+1) = N_extrinsic1 * W(j) + E2(t) * W_intrinsic - E2(t) * W(j) + (c + (b-c).*rand(1))  ; % E1
                rate1(t+1) = rate1(t) + sig(E1(t+1));

                E2(t+1) = N_extrinsic2 * W(j) + E1(t) * W_intrinsic - E1(t) * W(j) + (c + (b-c).*rand(1)) ; % E2
                rate2(t+1) = rate2(t) + sig(E2(t+1));       
                end       
            end % end of time
            if rate2(51) > rate1(51)
                tmp1(trial)= 1;
            end
            if rate2(102) > rate1(102)
                tmp2(trial)=1;
            end                             
                  
        end % end of trials 
        M50(i,j) = sum(tmp1); 
        M100(i,j) = sum(tmp2);
        
        if i==2 && j==3
        figure
        plot(time,rate1(2:102),'r')
        xlim([0 100])
        title(strcat('Activities of neurons E1 and E2 at W.extrinsic=',num2str(W(j))))
        xlabel('time')
        ylabel('firing rate')
        hold on 
        plot(time,rate2(2:102),'k')
        legend('E2','E1')
        hold off
        end
%         disp(j)
    end
% disp(i)
end

figure
subplot(3,2,1)
imagesc(W,a,M50,[0 100])
colorbar
title('Modulate Extrinsic and Inhibitory weights: Time Step 50')
xlabel('Extrinsic and Inhibitory Weights')
ylabel('Activity Difference')
set(gca,'XTickLabel',{'0.4','1.6','3.6','6.4','10'},'XTick',[0.4,1.6,3.6,6.4,10]);

subplot(3,2,2)
imagesc(W,a,M100,[0 100])
colorbar
title('Modulate Extrinsic and Inhibitory weights: Time Step 100')
xlabel('Extrinsic and Inhibitory Weights')
ylabel('Activity Difference')
xticklabels = [0.4,1.6,3.6,6.4,10];
xticks = [0.4,1.6,3.6,6.4,10];
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
hold on



a = [0.05:0.05:0.25]; % activity difference, calculated from EXT2-EXT1 
W = [0.4, 1.6, 3.6, 6.4, 10];
W_fixed = 0.1;
M50 = zeros(length(a),length(W));
M100 = zeros(length(a),length(W));

t=1; % reset t

for i = 1:length(a)
    
    for j = 1:length(W)
            tmp1= zeros(1,100);
            tmp2= zeros(1,100);
        for trial = 1:100

            
            for t = 2:length(time)
                N_extrinsic1=0.5; 
                N_extrinsic2 = a(i) + N_extrinsic1; 
                
                E1(t+1) = N_extrinsic1 * W_fixed + E2(t) * W_fixed - E2(t) * W_fixed + (c + (b-c).*rand(1))  ; % E1
                rate1(t+1) = rate1(t) + sig(E1(t+1));

                E2(t+1) = N_extrinsic2 * W_fixed + E1(t) * W_fixed - E1(t) * W_fixed + (c + (b-c).*rand(1)) ; % E2
                rate2(t+1) = rate2(t) + sig(E2(t+1));
                             
                if t > 52  % midtrial, 50th time step     
                E1(t+1) = N_extrinsic1 *W_fixed + E2(t) * W(j) - E2(t) *  + (c + (b-c).*rand(1))  ; % E1
                rate1(t+1) = rate1(t) + sig(E1(t+1));

                E2(t+1) = N_extrinsic2 *W_fixed + E1(t) * W(j) - E1(t) * W_fixed + (c + (b-c).*rand(1)) ; % E2
                rate2(t+1) = rate2(t) + sig(E2(t+1));       
                end       
            end % end of time
            if rate2(51) > rate1(51)
                tmp1(trial)= 1;
            end
            if rate2(102) > rate1(102)
                tmp2(trial)=1;
            end                             
                  
        end % end of trials 
        M50(i,j) = sum(tmp1); 
        M100(i,j) = sum(tmp2);

    end

end


subplot(3,2,3)
imagesc(W,a,M50,[0 100])
colorbar
title('Modulate Extrinsic and Inhibitory weights: Time Step 50')
xlabel('Extrinsic and Inhibitory Weights')
ylabel('Activity Difference')
set(gca,'XTickLabel',{'0.4','1.6','3.6','6.4','10'},'XTick',[0.4,1.6,3.6,6.4,10]);

subplot(3,2,4)
imagesc(W,a,M100,[0 100])
colorbar
title('Modulate Extrinsic and Inhibitory weights: Time Step 100')
xlabel('Extrinsic and Inhibitory Weights')
ylabel('Activity Difference')
xticklabels = [0.4,1.6,3.6,6.4,10];
xticks = [0.4,1.6,3.6,6.4,10];
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)









