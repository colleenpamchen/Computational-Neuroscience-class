function err = jk_hw1_popcoding (nPrefDirs)

% Generate population code using cosine tuning curves

% get a uniform distribution for the preferred direction of a cell and its
% maximum firing rate
for i = 1:nPrefDirs
    prefDir(i) = random('unif',0,2*pi);
    fr(i) = random('unif', 5, 50);
end

% test on eight different directions
inx = 0;
for i = pi/8:pi/8:2*pi
    inx = inx + 1;
    ang(inx) = i;
    
    % the error is the sine of the angle between the desired angle and the
    % angle calculated from the population code
    err(inx) = abs(sin(i - popCode(i,prefDir, fr)));
    
end

figure;
plot(ang, err);
axis ([0 2*pi 0 1])
xlabel ('Angle in Radians');
ylabel ('Error in Radians');

end

function popAng = popCode(a, pd, fr)

% a is the angle
% pd is an array of preferred directions
% fr is an array of maximum firing rates
% returns the angle of the population vector.

popX = 0;
popY = 0;

figure;

% show the desired angle in red using a compass plot
compass(cos(a),sin(a), 'r')
hold on

for i = 1:size(pd,2)
    
    % calculate the x and y components of the vector using the preferred
    % direction and the max firing rate of each population member. Plot the
    % member's contribution in cyan on the compass plot.
    x = cos(pd(i)) * max(0,fr(i)*cos(a-pd(i)));
    y = sin(pd(i)) * max(0,fr(i)*cos(a-pd(i)));
    compass(x/max(fr),y/max(fr), 'c')
    
    % Use vector addition to create a population code.
    popX = popX + x;
    popY = popY + y;
end

% Plot the vector population in blue on the compass plot
compass(popX/max(fr),popY/max(fr),'b')

% Use the arc tangent command to get the angle of the population vector.
popAng = atan2(popY, popX);

end



