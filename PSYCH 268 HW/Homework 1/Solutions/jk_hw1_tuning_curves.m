function jk_hw1_tuning_curves

% Part 1: Cosine Tuning Curves 

% change firing rate parameters
subplot(1,2,1)

% Replicate response in figure 1.6B
r0 = 32.34;     % baseline activity for 1.6B
rmax = 54.69;   % max firing rate
smax = 161.25;  % preferred direction in degrees
tw = 1;
cosineTuning (r0, rmax, smax, tw, 'k')
hold on

% Change the base rate
r0 = 15;     % half of example 1.6B
rmax = 54.69;   % max firing rate
smax = 160;  % preferred direction in degrees
tw = 1;
cosineTuning (r0, rmax, smax, tw, 'g')

% Change the max firing rate
r0 = 32.34;     % baseline activity for 1.6B
rmax = 100;   % max firing rate
smax = 162;  % preferred direction in degrees
tw = 1;
cosineTuning (r0, rmax, smax, tw, 'b')

legend ('1.6B', 'r0=15', 'rmax=100')
axis([0 360 0 110])
xlabel ('s(movement direction in degrees')
ylabel ('f(Hz)')

% change tuning widths
subplot(1,2,2)

% Replicate response in figure 1.6B
r0 = 32.34;     % baseline activity for 1.6B
rmax = 54.69;   % max firing rate
smax = 161.25;  % preferred direction in degrees
tw = 1;
cosineTuning (r0, rmax, smax, tw, 'k')
hold on

% Replicate response in figure 1.6B with no base rate
r0 = 0;     % baseline activity for 1.6B
rmax = 54.69;   % max firing rate
smax = 161.25;  % preferred direction in degrees
tw = 1;
cosineTuning (r0, rmax, smax, tw, 'r')
hold on

% Narrow the cosine tuning curve
r0 = 0;     % baseline activity for 1.6B
rmax = 54.69;   % max firing rate
smax = 161.25;  % preferred direction in degrees
tw = 5;
cosineTuning (r0, rmax, smax, tw, 'g')

% Even more narrow
r0 = 0;     % baseline activity for 1.6B
rmax = 54.69;   % max firing rate
smax = 161.25;  % preferred direction in degrees
tw = 9;
cosineTuning (r0, rmax, smax, tw, 'b')

legend ('1.6B', 'broad', 'narrow', 'narrower')
axis ([0 360 0 60])
xlabel ('s(movement direction in degrees')
ylabel ('f(Hz)')

% Part 2: Sigmoid Tuning Curve

figure;
subplot (1,2,1)

% replicate 1.7B
rmax = 36.03;
sHalf = 0.036;
deltaS = 0.029;
sigmoidTuning(rmax, sHalf, deltaS, 'k');

hold on

% Change the slope of the sigmoid
deltaS = 0.25;
sigmoidTuning(rmax, sHalf, deltaS, 'r');
deltaS = 0.005;
sigmoidTuning(rmax, sHalf, deltaS, 'g');

legend ('1.7B', 'ds=0.25', 'ds=0.005')
axis([-1 1 0 40])
xlabel ('s(retinal disparity in degrees')
ylabel ('f(Hz)')

subplot (1,2,2)

% replicate 1.7B
rmax = 36.03;
sHalf = 0.036;
deltaS = 0.029;
sigmoidTuning(rmax, sHalf, deltaS, 'k');

hold on

% change the max firing rate
rmax = 15;
sigmoidTuning(rmax, sHalf, deltaS, 'r');

% change the halfway point
rmax = 36.03;
sHalf = 0.5;
sigmoidTuning(rmax, sHalf, deltaS, 'g');

legend ('1.7B', 'rmax=15', 's1/2=0.5')
axis([-1 1 0 40])
xlabel ('s(retinal disparity in degrees')
ylabel ('f(Hz)')

end

function cosineTuning (r0, rmax, smax, tw, color)

%   convert tw to radians
smaxRad = smax/360*2*pi;

for i = 1:360
    sRad = i/360*2*pi;
    f(i) = max(0, r0 + (rmax - r0)*cos(sRad - smaxRad)^tw);
end

plot(1:360,f, color)

end

function sigmoidTuning (rmax, sHalf, deltaS, color)

s=-1:0.005:1;

for i = 1:size(s,2)
    f(i) = rmax/(1 + exp((sHalf - s(i))/deltaS));
end

plot (s,f, color)

end