function jk_hw3_phasicnm

defaultGAIN = 1;

T = 100;
TRIALS = 100;
ext1 = 0.50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('FIG4A - EXTRINSIC & INHIBITORY WGTS')
wgtPhasic = [0.1, 0.4, 0.8, 1.6, 2.5, 3.6, 4.8, 6.4, 8.2, 10.0];
actDiff = [0.05 0.10 0.15 0.20 0.25];

fig4a_50 = zeros(size(actDiff,2),size(wgtPhasic,2));
fig4a_100 = zeros(size(actDiff,2),size(wgtPhasic,2));

% loop through a range of activity differences and phasic weights
for a = 1:size(actDiff,2)
    for w = 1:size(wgtPhasic,2)
        disp(['   actDiff= ', num2str(actDiff(a)), ' wgtPhasic=', num2str(wgtPhasic(w))])
        ext2 = ext1 + actDiff(a);
        
        % run for many trials
        for t = 1:TRIALS
            e1 = zeros(1,T);
            e2 = zeros(1,T);
            sumE1 = zeros(1,T);
            sumE2 = zeros(1,T);
            wExt = 0.10;
            wInt = 0.10;
            wInh = 0.10;
            
            % each trial lasts 100 time steps
            for i = 2:T
                
                % get the synaptic input and run that through a sigmoid
                % function with a set gain.
                e1(i) = sigmoid(ext1*wExt + e2(i-1)*wInt - e2(i-1)*wInh + random('unif', -1, 1), defaultGAIN);
                e2(i) = sigmoid(ext2*wExt + e1(i-1)*wInt - e1(i-1)*wInh + random('unif', -1, 1), defaultGAIN);
                sumE1(i) = sumE1(i-1) + e1(i);
                sumE2(i) = sumE2(i-1) + e2(i);
                
                % at the halfway point simulate phasic neumodulation by
                % increasing the extrinsic and inhibitory weights
                if i == 50
                    wExt = wgtPhasic(w);
                    wInh = wgtPhasic(w);
                end
            end
            
            % calculate the cumulative firing rate
            if sumE2(50) > sumE1(50)
                fig4a_50(a,w) = fig4a_50(a,w) + 1;
            end
            if sumE2(100) > sumE1(100)
                fig4a_100(a,w) = fig4a_100(a,w) + 1;
            end
            
            % for the representative example, the extrinsic drive is 0.60
            % and pick an example where the final value of E2 is at least
            % twice that of E1
            if ext2 == 0.6 && sumE2(100) > 2*sumE1(100)
                plot(1:T,sumE1,1:T,sumE2)
                legend('E1', 'E2')
                xlabel('Time')
                ylabel('Cumulative Firing Rate')
            end
            
        end
    end
end

figure;

subplot(3,2,1)
imagesc(fig4a_50,[1 100])
axis square
colorbar
title('Modulation of Extrinsic and Inhibitory weights (timestep 50)')

subplot(3,2,2)
imagesc(fig4a_100,[1 100])
axis square
colorbar
title('Modulation of Extrinsic and Inhibitory weights (timestep 100)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('FIG4B - INTRINSIC WGTS')
wgtPhasic = [0.1, 0.4, 0.8, 1.6, 2.5, 3.6, 4.8, 6.4, 8.2, 10.0];
actDiff = [0.05 0.10 0.15 0.20 0.25];

fig4b_50 = zeros(size(actDiff,2),size(wgtPhasic,2));
fig4b_100 = zeros(size(actDiff,2),size(wgtPhasic,2));

% loop through a range of activity differences and phasic weights
for a = 1:size(actDiff,2)
    for w = 1:size(wgtPhasic,2)
        disp(['   actDiff= ', num2str(actDiff(a)), ' wgtPhasic=', num2str(wgtPhasic(w))])
        ext2 = ext1 + actDiff(a);
        
        % run for many trials       
        for t = 1:TRIALS
            e1 = zeros(1,T);
            e2 = zeros(1,T);
            sumE1 = zeros(1,T);
            sumE2 = zeros(1,T);
            wExt = 0.10;
            wInt = 0.10;
            wInh = 0.10;
            
            % each trial lasts 100 timesteps
            for i = 2:T
                % get the synaptic input and run that through a sigmoid
                % function with a set gain.
                e1(i) = sigmoid(ext1*wExt + e2(i-1)*wInt - e2(i-1)*wInh + random('unif', -1, 1), defaultGAIN);
                e2(i) = sigmoid(ext2*wExt + e1(i-1)*wInt - e1(i-1)*wInh + random('unif', -1, 1), defaultGAIN);
                sumE1(i) = sumE1(i-1) + e1(i);
                sumE2(i) = sumE2(i-1) + e2(i);
                
                % at the halfway point simulate phasic neuromodulation by
                % increasing the intrinsic weights
                if i == 50
                    wInt = wgtPhasic(w);
                end
            end
            
            % calculate the cumulative firing rate            
            if sumE2(50) > sumE1(50)
                fig4b_50(a,w) = fig4b_50(a,w) + 1;
            end
            if sumE2(100) > sumE1(100)
                fig4b_100(a,w) = fig4b_100(a,w) + 1;
            end
        end
    end
end

subplot(3,2,3)
imagesc(fig4b_50,[1 100])
axis square
colorbar
title('Modulation of Intrinsic (timestep 50)')

subplot(3,2,4)
imagesc(fig4b_100,[1 100])
axis square
colorbar
title('Modulation of Intrinsic (timestep 100)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('FIG4C - SYNAPTIC GAIN')
gainPhasic = [1 2 3 5 5 6 7 8 9 10];
actDiff = [0.05 0.10 0.15 0.20 0.25];

fig4c_50 = zeros(size(actDiff,2),size(wgtPhasic,2));
fig4c_100 = zeros(size(actDiff,2),size(wgtPhasic,2));

% loop through a range of activity differences and phasic synaptic gains
for a = 1:size(actDiff,2)
    for g = 1:size(wgtPhasic,2)
        disp(['   actDiff= ', num2str(actDiff(a)), ' wgtPhasic=', num2str(wgtPhasic(g))])
        ext2 = ext1 + actDiff(a);
        
        for t = 1:TRIALS
            e1 = zeros(1,T);
            e2 = zeros(1,T);
            sumE1 = zeros(1,T);
            sumE2 = zeros(1,T);
            wExt = 0.10;
            wInt = 0.10;
            wInh = 0.10;
            gain = defaultGAIN;
            
            for i = 2:T
                % get the synaptic input and run that through a sigmoid
                % function with a variable gain.
                e1(i) = sigmoid(ext1*wExt + e2(i-1)*wInt - e2(i-1)*wInh + random('unif', -1, 1), gain);
                e2(i) = sigmoid(ext2*wExt + e1(i-1)*wInt - e1(i-1)*wInh + random('unif', -1, 1), gain);
                sumE1(i) = sumE1(i-1) + e1(i);
                sumE2(i) = sumE2(i-1) + e2(i);
                
                % at the halfway point, simulate phasic neuromodulation by
                % increasing the synaptic gain.
                if i == 50
                    gain = gainPhasic(g);
                end
            end
            
            % calculate the cumulative firing rate
            if sumE2(50) > sumE1(50)
                fig4c_50(a,g) = fig4c_50(a,g) + 1;
            end
            if sumE2(100) > sumE1(100)
                fig4c_100(a,g) = fig4c_100(a,g) + 1;
            end
        end
    end
end

subplot(3,2,5)
imagesc(fig4c_50,[1 100])
axis square
colorbar
title('Modulation of Synaptic Gain (timestep 50)')

subplot(3,2,6)
imagesc(fig4c_100,[1 100])
axis square
colorbar
title('Modulation of Synaptic Gain (timestep 100)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function s = sigmoid(synIn, gain)

s = 1/(1 + exp(-gain*synIn));

end