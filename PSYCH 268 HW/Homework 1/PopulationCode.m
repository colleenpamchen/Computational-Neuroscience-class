% Population Neuron Coding
% Generate a population code of direction using cosine tuning curves(see section 3.3). 
% – Vary the population size(e.g.,4,8,16,32) and preferred directions to see how this affects the coding error.

figure;
hold on
rng('shuffle')

%% Population Size 4
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

%% Population Size 8
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

%% Population Size 16
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

%% Population Size 32
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

%% Plots
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

