%% Izhikevich Parameter 1
clear;
clc;
figure;
hold on;

a = 0.02;
b = -0.1;
c = -55;
d = 6;

vr = 0; % arbitrary resting potential
vt = 30; % arbitrary threshold potential

for I = 1:30
    v = c;
    u = b * c;
    count = 0;
    
    for k = 1:50
        vi(I,k) = (k * ((v-c)^2)) - u + I;
        ui(I,k) = a * ((b * v) - u);
        
        v = v + vi(I,k);
        u = u + ui(I,k);
        
        if vi(I,k) > vt
            v = c;
            u = u + d;
            count = count + 1; 
        end
        vstore(I,k) = v;
        counter1(I) = count;
    end 
end

%% Izhikevich Parameter 2

a = 0.2;
b = 0.26;
c = -65;
d = 0;

vr = 0; % arbitrary resting potential
vt = 30; % arbitrary threshold potential

for I = 1:30
    v = c;
    u = b * c;
    count = 0;
    
    for k = 1:50
        vi(I,k) = (k * ((v-c)^2)) - u + I;
        ui(I,k) = a * ((b * v) - u);
        
        v = v + vi(I,k);
        u = u + ui(I,k);
        
        if vi(I,k) > vt
            v = c;
            u = u + d;
            count = count + 1; 
        end
        vstore(I,k) = v;
        counter2(I) = count;
    end 
end

%% Izhikevich Parameter 3

a = 0.02;
b = 0.2;
c = -65;
d = 8;

vr = 0; % arbitrary resting potential
vt = 30; % arbitrary threshold potential

for I = 1:30
    v = c;
    u = b * c;
    count = 0;
    
    for k = 1:1000
        vi(I,k) = (0.04 * ((v-c)^2)) - u + I;
        ui(I,k) = a * ((b * v) - u);
        
        v = v + vi(I,k);
        u = u + ui(I,k);
        
        if v > vt
            v = c;
            u = u + d;
            count = count + 1; 
        end
        vstore(I,k) = v;
        counter3(I) = count;
    end 
end

plot(1:30, counter1, 1:30, counter2, 1:30, counter3)
legend('Param 1', 'Param 2', 'Param 3')
xlabel('Current')
ylabel('Frequency')
title('Izhikevich Neurons')