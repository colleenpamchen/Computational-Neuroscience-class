% Veronica Chu
% Psych 268 Homework 1 - Oct 6, 2015
% Note: mean error is incorrect, but this was as close as I was
% able to get it
%% Cosine Tuning Curve (M1 neuron)
clear;
clc;
figure;
hold on

ro = 32.34;
rmax = 54.69;
smax = 161.25;

% Original Tuning Curve
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

% Narrower Tuning Curve
avgresp2 = NaN(36,1);
s = 10;
while s < 40
    for i = 1:36
        avgresp2(i) = ro + (rmax - ro)*cos(2*(s - smax));
        if avgresp2(i) < 0
            avgresp2(i) = 0;
        end
    s = s + 0.2;
    end
end

% Wider Tuning Curve
avgresp3 = NaN(36,1);
s = 10;
while s < 40
    for i = 1:36
        avgresp3(i) = ro + (rmax - ro)*cos(0.5*(s - smax));
        if avgresp3(i) < 0
            avgresp3(i) = 0;
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

subplot(1,3,2)
plot(avgresp2)
xlim([0 36])
ylim([0 60])
set(gca,'XTickLabel',[0 50 100 150 200 250 300 360]);
title('M1 Neuron - Narrower Curve')
xlabel('s (movement direction in degrees)')
ylabel('f (Hz)')

subplot(1,3,3)
plot(avgresp3)
xlim([0 36])
ylim([0 60])
set(gca,'XTickLabel',[0 50 100 150 200 250 300 360]);
title('M1 Neuron - Wider Curve')
xlabel('s (movement direction in degrees)')
ylabel('f (Hz)')

%% Sigmoidal Tuning Curve (V1 neuron with retinal disparity)
clear;
clc;
figure;
hold on

rmax = 36.03;
shalf = 0.036;
acc = 0.029;

% Original Tuning Curve
avgresp1 = NaN(10,1);
s = -1.0;
while s < 1.5
    for i = 1:10
        avgresp1(i) = rmax / (1 + exp((shalf - s)/acc));
        s = s + 0.3;
    end
end

% Lower Gain Tuning Curve
avgresp2 = NaN(10,1);
s = -1.0;
while s < 1.5
    for i = 1:10
        avgresp2(i) = rmax / (2*(1 + exp((shalf - s)/acc)));
        s = s + 0.3;
    end
end

% Higher Gain Tuning Curve
avgresp3 = NaN(10,1);
s = -1.0;
while s < 1.5
    for i = 1:10
        avgresp3(i) = (2 * rmax) / (1 + exp((shalf - s)/acc));
        s = s + 0.3;
    end
end

subplot(1,3,1)
plot(avgresp1)
xlim([1 7])
set(gca,'XTickLabel',[-1.5 -1.0 -0.5 0 0.5 1.0 1.5]);
title('V1 Neuron (with Retinal Disparity) - Original Curve')
xlabel('s (retinal disparity in degrees)')
ylabel('f (Hz)')

subplot(1,3,2)
plot(avgresp2)
xlim([1 7])
set(gca,'XTickLabel',[-1.5 -1.0 -0.5 0 0.5 1.0 1.5]);
title('V1 Neuron (with Retinal Disparity) - Lower Gain Curve')
xlabel('s (retinal disparity in degrees)')
ylabel('f (Hz)')

subplot(1,3,3)
plot(avgresp3)
xlim([1 7])
set(gca,'XTickLabel',[-1.5 -1.0 -0.5 0 0.5 1.0 1.5]);
title('V1 Neuron (with Retinal Disparity) - Higher Gain Curve')
xlabel('s (retinal disparity in degrees)')
ylabel('f (Hz)')

%% Population Neuron Coding
clear;
clc;
figure;
hold on
rng('shuffle')

% Population Size 4
c = randi(360,1,4);
rmax = randi(100,1,4);
pop4 = NaN(90,4);

for a = 1:4
    for r = 1:360
        pop4(r,a) = (rmax(a))*cos(0.015*(r - c(a)));
        if pop4(r,a) < 0
           pop4(r,a) = 0; 
        end
    end
end

wind = randi(360,1,3);  % actual wind direction
popang4 = NaN(3,4);
r = NaN(3,4);
popvec4 = NaN(1,3);
meanerror4 = NaN(1,3);

for i = 1:3
    for a = 1:4
        r(i,a) = (wind(i) * rmax(a)) / c(a); % response from actual input
        popang4(i,a) = (r(i,a)/rmax(a)) * c(a); % population vector
    end
end
for i = 1:3
    popvec4(i) = sum(popang4(i,:));
    meanerror4(i) = sqrt(popvec4(i).^2 - r(i).^2);
    meanerror4(i) = meanerror4(i)/360;
end

% Population Size 8
c = randi(360,1,8);
rmax = randi(100,1,8);
pop8 = NaN(90,8);

for a = 1:8
    for r = 1:360
        pop8(r,a) = (rmax(a))*cos(0.015*(r - c(a)));
        if pop8(r,a) < 0
           pop8(r,a) = 0; 
        end
    end
end

wind = randi(360,1,3);  % actual wind direction
popang8 = NaN(3,8);
r = NaN(3,4);
popvec8 = NaN(1,3);
meanerror8 = NaN(1,3);

for i = 1:3
    for a = 1:8
        r(i,a) = (wind(i) * rmax(a)) / c(a); % response from actual input
        popang8(i,a) = (r(i,a)/rmax(a)) * c(a); % population vector
    end
end
for i = 1:3
    popvec8(i) = sum(popang8(i,:));
    meanerror8(i) = sqrt(popvec8(i).^2 - r(i).^2);
    meanerror8(i) = meanerror8(i)/360;
end

% Population Size 16
c = randi(360,1,16);
rmax = randi(100,1,16);
pop16 = NaN(90,16);

for a = 1:16
    for r = 1:360
        pop16(r,a) = (rmax(a))*cos(0.015*(r - c(a)));
        if pop16(r,a) < 0
           pop16(r,a) = 0; 
        end
    end
end

wind = randi(360,1,3);  % actual wind direction
popang16 = NaN(3,16);
r = NaN(3,16);
popvec16 = NaN(1,3);
meanerror16 = NaN(1,3);

for i = 1:3
    for a = 1:16
        r(i,a) = (wind(i) * rmax(a)) / c(a); % response from actual input
        popang16(i,a) = (r(i,a)/rmax(a)) * c(a); % population vector
    end
end
for i = 1:3
    popvec16(i) = sum(popang16(i,:));
    meanerror16(i) = sqrt(popvec16(i).^2 - r(i).^2);
    meanerror16(i) = meanerror16(i)/360;
end

% Population Size 32
c = randi(360,1,32);
rmax = randi(100,1,32);
pop32 = NaN(90,32);

for a = 1:32
    for r = 1:360
        pop32(r,a) = (rmax(a))*cos(0.015*(r - c(a)));
        if pop32(r,a) < 0
           pop32(r,a) = 0; 
        end
    end
end

wind = randi(360,1,3);  % actual wind direction
popang32 = NaN(3,32);
r = NaN(3,32);
popvec32 = NaN(1,3);
meanerror32 = NaN(1,3);

for i = 1:3
    for a = 1:32
        r(i,a) = (wind(i) * rmax(a)) / c(a); % response from actual input
        popang32(i,a) = (r(i,a)/rmax(a)) * c(a); % population vector
    end
end
for i = 1:3
    popvec32(i) = sum(popang32(i,:));
    meanerror32(i) = sqrt(popvec32(i).^2 - r(i).^2);
    meanerror32(i) = meanerror32(i)/360;
end

% Plots
subplot(2,4,1), plot(pop4)
title('Population Code - 4 Neurons')
subplot(2,4,2), plot(pop8)
title('Population Code - 8 Neurons')
subplot(2,4,3), plot(pop16)
title('Population Code - 16 Neurons')
subplot(2,4,4), plot(pop32)
title('Population Code - 32 Neurons')

subplot(2,4,5)
plot(wind,meanerror4,'o')
xlim([0 360])
title('Mean Error - 4 Neurons')
xlabel('Wind Direction (degrees)')
ylabel('Error (degrees)')

subplot(2,4,6)
plot(wind,meanerror8,'o')
xlim([0 360])
title('Mean Error - 8 Neurons')
xlabel('Wind Direction (degrees)')
ylabel('Error (degrees)')

subplot(2,4,7)
plot(wind,meanerror16,'o')
xlim([0 360])
title('Mean Error - 16 Neurons')
xlabel('Wind Direction (degrees)')
ylabel('Error (degrees)')

subplot(2,4,8)
plot(wind,meanerror32,'o')
xlim([0 360])
title('Mean Error - 32 Neurons')
xlabel('Wind Direction (degrees)')
ylabel('Error (degrees)')

