% Veronica Chu 
% Psych 268 Homework 3 - Oct 19, 2015

figure;
hold on

%% Preferred Orientations
c = [45 90 135 180];
rmax = [50 50 50 50];
pop4 = NaN(90,4);

for a = 1:4
    for r = 1:360
        pop4(r,a) = (rmax(a))*cos(0.015*(r - c(a)));
        if pop4(r,a) < 0
           pop4(r,a) = 0; 
        end
    end
end

subplot(1,3,1), plot(pop4)
title('Preferred Orientations of Neurons')
xlabel('Orientation (deg)')
ylabel('Firing Rate (Hz)')
xlim([0 360])

%% Horizontal and Vertical

w = [1 4 1 2];  % adjust weights for each orientation
c = [45 90 135 180];
rmax = [50 50 50 50];
pop4 = NaN(90,4);

for r = 1:360
    pop4(r,1) = w(1) * 0;
    
    pop4(r,2) = w(2) * (rmax(2))*cos(0.015*(r - c(2)));
    if pop4(r,2) < 0
        pop4(r,2) = 0;
    end
    
    pop4(r,3) = w(3) * 0;
    
    pop4(r,4) = w(4) * (rmax(4))*cos(0.015*(r - c(4)));
    if pop4(r,4) < 0
        pop4(r,4) = 0;
    end
end

subplot(1,3,2), plot(pop4)
title('Vertical and Horizontal (No Inhibition)')
xlabel('Orientation (deg)')
ylabel('Firing Rate (Hz)')
xlim([0 360])

%% Horizontal and Vertical with Inhibition

w = [1 4 1 2];  % adjust weights for each neuron
c = [45 90 135 180];
rmax = [50 50 50 50];
pop4 = NaN(90,4);

for r = 1:360
    pop4(r,1) = w(1) * 0;
    pop4(r,2) = w(2) * (rmax(2))*cos(0.015*(r - c(2)));
    pop4(r,3) = w(3) * 0;
    pop4(r,4) = w(4) * (rmax(4))*cos(0.015*(r - c(4)));
    
    pop4(r,2) = max(pop4(r,2),0);
    pop4(r,4) = max(pop4(r,4),0);
     
    if pop4(r,2) > 0 && pop4(r,4) > 0 && pop4(r,2) > pop4(r,4)
        pop4(r,2) = pop4(r,2) - pop4(r,4);
        pop4(r,4) = 0;
    end
    
    if pop4(r,4) > 0 && pop4(r,2) > 0 && pop4(r,4) > pop4(r,2)
        pop4(r,4) = pop4(r,4) - pop4(r,2);
        pop4(r,2) = 0;
    end
    
    pop4(r,2) = max(pop4(r,2),0);
    pop4(r,4) = max(pop4(r,4),0);
end

subplot(1,3,3), plot(pop4)
title('Vertical and Horizontal (with Inhibition)')
xlabel('Orientation (deg)')
ylabel('Firing Rate (Hz)')
xlim([0 360])