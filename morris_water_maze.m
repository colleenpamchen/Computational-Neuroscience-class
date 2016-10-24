function [latency, pc, w, z] = morris_water_maze

global radius;
% global obstacle;

dirs = 8; % possible headings
beta = 2; % for softmax
radius = 1.0; % arena is 2 meters wide
sigma = 0.16; % place cell tuning width of 0.16m
reward_value = 1;
eta = 0.01; % learning rate

% set up place cell indices across circular arena
inx = 0;
for i = -radius:sigma/2:radius
    for j = -radius:sigma/2:radius
        if norm([i j]) < radius
            inx = inx + 1;
            pc.x(inx) = i;
            pc.y(inx) = j;
        end
    end
end

reward = [0.50 -0.50];
% obstacle = [1.0 1.0];

z = zeros(dirs,inx); % actor weights
w = zeros(1,inx); % critic weights
TRIALS = 32;
latency = zeros(1,TRIALS);

for trial = 1:TRIALS
    
    % get rat's initial position. start each trial in a different quadrant
    if mod(trial,4) == 1
        rat.x = -0.9;
        rat.y = 0;
    elseif mod(trial,4) == 2
        rat.x = 0;
        rat.y = 0.9;
    elseif mod(trial,4) == 3
        rat.x = 0.9;
        rat.y = 0;
    else
        rat.x = 0;
        rat.y = -0.9;
    end
    
    % let the rat explore for 100 time steps or until it gets reward
    t = 1;
    found_reward = 0;
    v = 0;
    a = zeros(1,8);
    
    % run for 100 moves or until a reward is found. whichever comes first
    while t <= 250 && ~found_reward;
        
        % choose an action and move rat to new location
        act = action_select (a, beta);
        rat = move (act, rat);
        
        % save off the old critic value. calculate the new critic value
        % based on the rat's current position. calculate the actor value
        % for each direction.
        vpre = v;
        v = 0;
        a = zeros(1,8);
        for i=1:inx
            
            % get place cell activity
            pc.act(i) = place_cell([rat.x rat.y], [pc.x(i) pc.y(i)]);
            
            % calculate current value for critic based on weight and place
            % cell activity
            v = v + w(i)*pc.act(i);
            
            % calculate actor values based on weights and place cell
            % activity
            for j = 1:dirs
                a(j) = a(j) + z(j,i)*pc.act(i);
            end
        end
        
        % if rat is within 0.2 meters of platform, give reward
        found_reward = norm([rat.x rat.y]-reward) < 0.2;
        %         if found_reward
        %             disp(['Found reward at time ', num2str(t), ' on trial ', num2str(trial)])
        %         end
        
        % get TD delta rule
        d = delta(v, vpre, reward_value * found_reward);
        
        % update the critic's weights
        w = w + eta*d*pc.act;
        
        % update the actor's weights
        z(act,:) = z(act,:) + eta*d*pc.act;
        
        t = t + 1;
        
        subplot(131)
        scatter3(pc.x,pc.y,pc.act)
        axis square;
        subplot(132)
        scatter(pc.x,pc.y,10,'b.')
        hold on
        scatter(rat.x, rat.y,50,'ro')
        scatter(reward(1),reward(2),100,'go')
        %	scatter(obstacle(1),obstacle(2),300,'co')
        if found_reward
            title(['Found reward at time ', num2str(t), ' on trial ', num2str(trial), '!!!'])
        else
            title(['Trial ', num2str(trial), ' at time ', num2str(t)])
        end
        axis square;
        hold off
        subplot(133)
        vectorPlot(z,pc);
        axis square;
        drawnow;
    end
    if found_reward
        disp(['Found reward at time ', num2str(t), ' on trial ', num2str(trial)])
        pause(1)
    else
        disp(['Timed out on trial ', num2str(trial)])
    end
    latency(trial) = t;
end
end

% place cell response is a 2D Gaussian based on the distance from the agent
% to the place cell center
function y = place_cell (loc, ctr)

w = 2*0.25^2;
y = exp(-norm(loc-ctr)^2/w);
end

% delta rule for temporal difference learning
function y = delta(vnext, v, rwd)

y = rwd + vnext - v;
end

% choose an action based on the actor weights
function act = action_select (m, beta)

% get softmax probabilities
for i = 1:size(m,2)
    p(i) = exp(beta*m(i))/sum(exp(beta*m));
end

%
r = rand; % random number for action selection
sumprob = 0;

% choose an action based on the probability distribution
i = 1;
done = 0;
while ~done && i <= size(m,2)
    sumprob = sumprob + p(i);
    if sumprob > r
        act = i;
        done = 1;
    end
    i = i + 1;
end

% disp ([r sumprob act p m])
end

function p = move(action, pos)

global radius;
global obstacle
speed = 0.1;

% get the heading of the agent based on the desired action
switch action
    case 1
        dir = 0;
    case 2
        dir = pi/4;
    case 3
        dir = pi/2;
    case 4
        dir = 3*pi/4;
    case 5
        dir = pi;
    case 6
        dir = 5*pi/4;
    case 7
        dir = 3*pi/2;
    case 8
        dir = 7*pi/4;
end

% get a new position for the agent
tmp.x = pos.x + speed*cos(dir);
tmp.y = pos.y + speed*sin(dir);

% check that the agent is within the boundaries of the arena. If it is not
% randomly shift the heading slightly until it is within bounds
while norm([tmp.x tmp.y]) > radius || norm([tmp.x tmp.y]-obstacle) < 0.5
    dir = dir + sign(rand-0.5)*pi/8;
    tmp.x = pos.x + speed*cos(dir);
    tmp.y = pos.y + speed*sin(dir);
end

p = tmp;
end

function vectorPlot (wgts, pc)

dirs = size(wgts,1);
inx = 1;
dx = zeros(sqrt(size(wgts,2)));
dy = zeros(sqrt(size(wgts,2)));

for i = 1:sqrt(size(wgts,2))
    for j = 1:sqrt(size(wgts,2))
        for k = 1:dirs
            dx(i,j) = dx(i,j) + wgts(k,inx) * cos(2*pi*(k/dirs));
            dy(i,j) = dy(i,j) + wgts(k,inx) * sin(2*pi*(k/dirs));
            x(i,j) = pc.x(inx);
            y(i,j) = pc.y(inx);
        end
        inx = inx + 1;
    end
end

% [x,y] = meshgrid(1:size(dx,2));
% y = y .*2 ./22 -1;
% x = x .*2 ./22 -1;
quiver(x,y,dx,dy)
axis([-1 1 -1 1])
end
























