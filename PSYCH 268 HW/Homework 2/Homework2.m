% Veronica Chu
% Psych 268 Homework 1 - Oct 13, 2015

%% Sodium Permeability and Membrane Potential
figure;
ko = 5;
nao = 120;
clo = 125;

ki = 12;
nai = 125;
cli = 5;
c = 1;

Em = NaN(1,100);
for i=1:100
    b = 0.02;
    if 44 < i && i < 56
        b = 20;
    end
    ions = (ko + b*nao + c*cli) / (ki + b*nai + c*clo);
    Em(i) = 58*log(ions);
end

plot(1:100, Em)
xlabel('Time')
ylabel('Membrane Potential')
title('Permeability and Potential')

%% Integrate and Fire Neuron
clear;
clc;
figure;

rest = -65;
Vth = -50;
V = rest;
timeconst = 10;
resist = 10;

delta = NaN(4,100);
Vstore = NaN(4,100);
counter = NaN(1,4);
for Ie = 1:4
    count = 0;
    for i = 1:100
        delta(Ie,i) = (rest - V + (resist * Ie)) / timeconst;
        V = V + delta(Ie,i);
        if V >= Vth
           V = -65;
           count = count + 1; 
        end
        Vstore(Ie,i) = V;
    end
    counter(Ie) = count;
end

plot(1:4,counter)
xlabel('Current')
ylabel('Frequency of Spikes')
title('Integrate and Fire Neuron')

%% Izhikevich Neuron

% Parameter 1
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

% Parameter 2
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

% Parameter 3
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