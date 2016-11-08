clear;
clc;

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