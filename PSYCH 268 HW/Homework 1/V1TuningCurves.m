% Sigmoidal Tuning Curve (V1 neuron with retinal disparity)
% Reproduce the tuning curves in 1.7B. 
% – Vary the gain(1.7B) of the curves.

%% Original Tuning Curve

clear;
clc;
figure;
hold on

rmax = 36.03;
shalf = 0.036;
acc = 0.029;

avgresp1 = NaN(10,1);
s = -1.0;
while s < 1.5
    for i = 1:10
        avgresp1(i) = rmax / (1 + exp((shalf - s)/acc));
        s = s + 0.3;
    end
end

subplot(1,3,1)
plot(avgresp1)
xlim([1 7])
ylim([0 40])
set(gca,'XTickLabel',[-1.5 -1.0 -0.5 0 0.5 1.0 1.5]);
title('V1 Neuron (with Retinal Disparity) - Original Curve')
xlabel('s (retinal disparity in degrees)')
ylabel('f (Hz)')

%% Lower Gain Tuning Curve

rmax = 36.03;
shalf = 0.036;
acc = 0.029;

avgresp2 = NaN(10,1);
s = -1.0;
while s < 1.5
    for i = 1:10
        avgresp2(i) = rmax / (2*(1 + exp((shalf - s)/acc)));
        s = s + 0.3;
    end
end

subplot(1,3,2)
plot(avgresp2)
xlim([1 7])
% ylim([0 40])
set(gca,'XTickLabel',[-1.5 -1.0 -0.5 0 0.5 1.0 1.5]);
title('V1 Neuron (with Retinal Disparity) - Lower Gain Curve')
xlabel('s (retinal disparity in degrees)')
ylabel('f (Hz)')


%% Higher Gain Tuning Curve

rmax = 36.03;
shalf = 0.036;
acc = 0.029;

avgresp3 = NaN(10,1);
s = -1.0;
while s < 1.5
    for i = 1:10
        avgresp3(i) = (2 * rmax) / (1 + exp((shalf - s)/acc));
        s = s + 0.3;
    end
end

subplot(1,3,3)
plot(avgresp3)
xlim([1 7])
% ylim([0 40])
set(gca,'XTickLabel',[-1.5 -1.0 -0.5 0 0.5 1.0 1.5]);
title('V1 Neuron (with Retinal Disparity) - Higher Gain Curve')
xlabel('s (retinal disparity in degrees)')
ylabel('f (Hz)')