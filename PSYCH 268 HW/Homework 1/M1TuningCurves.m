% Cosine Tuning Curve (M1 neuron)
% Reproduce the tuning curves in 1.6B. 
% – Vary the tuning width(1.6B) of the curves.
clear;
clc;
figure;
hold on

ro = 32.34;
rmax = 54.69;
smax = 161.25;

%% Original Tuning Curve
avgresp1 = NaN(36,1);
s = 10;
while s < 40
    for i = 1:36
        avgresp1(i) = ro + (rmax - ro)*cos(s - smax);
        if avgresp1(i) < 0
            avgresp1(i) = 0;
        end
        s = s + 0.2;
    end
end

subplot(1,3,1)
plot(avgresp1)
xlim([0 36])
ylim([0 60])
set(gca,'XTickLabel',[0 50 100 150 200 250 300 360]);
title('M1 Neuron - Original Curve')
xlabel('s (movement direction in degrees)')
ylabel('f (Hz)')

%% Narrower Tuning Curve
ro = 32.34;
rmax = 54.69;
smax = 161.25;

avgresp2 = NaN(36,1);
% s = 66;
s = 10;
% while s < 68
while s < 40
    for i = 1:36
        avgresp2(i) = ro + (rmax - ro)*cos(2*(s - smax));
        if avgresp2(i) < 0
            avgresp2(i) = 0;
        end
%         s = s + 0.0056;
    s = s + 0.2;
    end
end

subplot(1,3,2)
plot(avgresp2)
xlim([0 36])
ylim([0 60])
set(gca,'XTickLabel',[0 50 100 150 200 250 300 360]);
title('M1 Neuron - Narrower Curve')
xlabel('s (movement direction in degrees)')
ylabel('f (Hz)')

%% Wider Tuning Curve
ro = 32.34;
rmax = 54.69;
smax = 161.25;

avgresp3 = NaN(36,1);
% s = 20;
s = 10;
% while s < 35
while s < 40
    for i = 1:36
        avgresp3(i) = ro + (rmax - ro)*cos(0.5*(s - smax));
        if avgresp3(i) < 0
            avgresp3(i) = 0;
        end
        s = s + 0.2;
%         s = s + 0.03;
    end
end

subplot(1,3,3)
plot(avgresp3)
xlim([0 36])
ylim([0 60])
set(gca,'XTickLabel',[0 50 100 150 200 250 300 360]);
title('M1 Neuron - Wider Curve')
xlabel('s (movement direction in degrees)')
ylabel('f (Hz)')
