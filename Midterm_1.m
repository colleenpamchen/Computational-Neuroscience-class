%% Midterm Assignment

% Make gradient 
freq = 8;
angle = 90; 
Period= 128;
xx = freq * cosd(angle);
yy = freq * sind(angle);
[X Y] = meshgrid(1:Period);
grating = cos(2*pi*((xx*X/Period) +(yy*Y/Period)));
imshow(grating); 

C = 50; % contrast 
0.5* C*X


