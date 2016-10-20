%
% Copyright (c) 2011 Regents of the University of California. All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%
% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%
% 3. The names of its contributors may not be used to endorse or promote
%    products derived from this software without specific prior written
%    permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% CARL Roomba Network
% 
% Accompanies the paper "A Neurorobotic Platform to Test the Influence of
% Neuromodulatory Signaling on Anxious and Curious Behavior", submitted to
% Frontiers in Neurorobotics.
%
% Jeff Krichmar
% University of California, Irvine
% November 2012
function [logrnet, frames] = krichmar_frontiers2013_carlroomba (tc_5ht, tc_da, mpfc, ofc, s)

% behavior states
global STATE_WALL_FOLLOW;
global STATE_OPEN_FIELD;
global STATE_EXPLORE_OBJECT;
global STATE_FIND_HOME;
global STATE_AT_HOME;
global STATE_LEAVE_HOME;
STATE_OPEN_FIELD = 1;
STATE_EXPLORE_OBJECT = 2;
STATE_WALL_FOLLOW = 3;
STATE_FIND_HOME = 4;
STATE_AT_HOME = 5;
STATE_LEAVE_HOME = 6;

% events
global EVENT_OBJECT;
global EVENT_LIGHT;
global EVENT_BUMP;
global E;
E = 3;
EVENT_OBJECT =   1;
EVENT_LIGHT = 2;
EVENT_BUMP =    3;
event = zeros (1,E);

% neural network parameters
global N;
global NM;
global n;
global nm;
global achne;
global tonic;


% parameters for return home
FIND_HOME_TIMEOUT = 10;
global dock_time;
dock_time = 0;
global found_home;
found_home = 0;

% parameters for wall following
global WALL_LEFT;
global WALL_RIGHT;
global side;
WALL_LEFT = -1;
WALL_RIGHT = 1;
side = WALL_RIGHT;

% threshold for camera detecting bright light
LIGHT_THRESHOLD = 0.50;

% parameters for bump event
TOO_CLOSE = 200; % millimeters

% parameters for beam event
force_field = 242;
buoy_and_force_field = 250;

rand('seed',s)

% initialize Roomba
serRoombaPort=RoombaInit ('/dev/ttyUSB0');

% initialize laser range finder
serURGPort = serial('/dev/ttyACM0','BaudRate', 115200);
fopen (serURGPort);
ranges = urg_getdata(serURGPort);   % get an initial reading

% take some initial sensor readings.
for i=1:3
    AngleSensorRoomba(serRoombaPort);
    DistanceSensorRoomba(serRoombaPort);
end

% initialize neural network
roomba_net_init (tc_5ht, tc_da, mpfc, ofc);

vid = videoinput('linuxvideo',1,'YUYV_160x120')
vid.ReturnedColorSpace = 'grayscale';
preview(vid);

tic;    % start timer

behave_state_time = toc;
behave_state = STATE_WALL_FOLLOW;
new_behave_state = behave_state;
start_demo = 1;


current_time = toc;
loginx = 0;
lite = 0;

% run the network for approximately three minutes
for t = 1:240
    
    current_time = toc;
    
    % get urg laser information
    oldranges = ranges; % save last scan
    ranges = urg_getdata(serURGPort);
    rnglen = size(ranges,2);
    objdist = urg_findobj (ranges);
    
    % find the min, max and average distances to the right, middle, and
    % left of the roomba
    [right_min, right_avg, right_max] = urg_getdist(ranges(round(rnglen*2/8)+1:round(rnglen*3/8)));
    [middle_min, middle_avg, middle_max] = urg_getdist(ranges(round(rnglen*3/8)+1:round(rnglen*5/8)));
    [left_min, left_avg, left_max] = urg_getdist(ranges(round(rnglen*5/8)+1:round(rnglen*6/8)));
    
    % get bump sensor information
    [BumpRight,BumpLeft,WheDropRight,WheDropLeft,WheDropCaster,BumpFront] = BumpsWheelDropsSensorsRoomba(serRoombaPort);
    
    % get dock beam status
    beam = DockBeam(serRoombaPort);
    homed = beam == force_field || beam >= buoy_and_force_field; % homed if charging or close to the docking station
    
    % grab a frame from the camera
    camframe = getsnapshot(vid);
    
    % get events (binary 1 == event occurred)
    lite = 0.0*lite + 1.0*(sum(sum(camframe,1),2)/(120*160*256)); % smooth the light value
    event(EVENT_LIGHT) = lite > LIGHT_THRESHOLD;
    if event(EVENT_LIGHT)
        disp([num2str(t), ': Light Event']);
    end
    event(EVENT_BUMP) = BumpLeft || BumpRight || BumpFront || min(min(left_min,right_min),middle_min) < TOO_CLOSE;
    event(EVENT_OBJECT) = objdist > 0; % object found in laser scan
    
    [chc] = roomba_net_cycle (event);
    if chc
        new_behave_state = chc;
    end
    
    % transitioned to a new state, print state information
    if behave_state ~= new_behave_state
        behave_state_time = toc;
        fwrite(serRoombaPort, [132]);   % go into full mode
        start_demo = 1;
        
        % set the side to turn when bumped (used by multiple states)
        if rand < 0.5
            side = WALL_LEFT;
        else
            side = WALL_RIGHT;
        end
        
        switch new_behave_state
            case STATE_WALL_FOLLOW
                disp ('WallFollow');
            case STATE_OPEN_FIELD
                disp ('OpenField');
            case STATE_EXPLORE_OBJECT
                disp ('ExploreObject');
            case STATE_FIND_HOME
                disp ('FindHome');
            case STATE_AT_HOME
                dock_time = toc;
                disp ('AtHome');
            case STATE_LEAVE_HOME
                disp ('LeaveHome');
        end
        behave_state = new_behave_state;
    end
    
    % handle states
    switch behave_state
        case STATE_WALL_FOLLOW
            WallFollow (serRoombaPort, start_demo);
        case STATE_OPEN_FIELD
            OpenField (serRoombaPort, event(EVENT_BUMP), left_avg, middle_avg, right_avg, 0.20, 0.40);
        case STATE_EXPLORE_OBJECT
            ExploreObject (serRoombaPort, event(EVENT_BUMP), objdist, rnglen, left_avg, middle_avg, right_avg);
        case STATE_FIND_HOME
            FindHome (serRoombaPort, start_demo);
        case STATE_AT_HOME
            AtHome (serRoombaPort);
        case STATE_LEAVE_HOME
            LeaveHome (serRoombaPort);
    end
    start_demo = 0;
    
    % logrnet state, event, and neural network information for post-processing
    loginx = loginx + 1;
    logrnet(loginx, 1) = current_time;
    logrnet(loginx, 2) = behave_state;
    logrnet(loginx, 3:3+E-1) = event;
    logrnet(loginx, 3+E:3+E+N+NM+E-1) = [n nm achne];
    logrnet(loginx, 3+E:3+E+N+NM+E+NM-1) = [n nm achne tonic];
    frames(loginx,:,:,:) = camframe;
    
end

% trial completed. stop roomba and close ports
fwrite(serRoombaPort, [132]);   % go into full mode
SetDriveWheelsCreate(serRoombaPort, 0.0, 0.0);
pause(1);
fclose (serURGPort);
fclose (serRoombaPort);
delete(vid);
save(['log-', num2str(tc_5ht), '-', num2str(tc_da), '-', num2str(mpfc), '-', num2str(ofc), '-', num2str(s), '.mat'], 'logrnet', 'frames')

end

% WallFollow (behavior state)
%
% Description
%   Attempts to keep roomba near a wall or large object based on laser
%   readings. If a collision is detected, it turns away from the collision.
%
% Inputs
%   serPort - serial port file id for roomba
%   bump - 1 if collision detected
%   dist - distance in mm detected by the laser range finder.
%
% Outputs
%   none
%
% function WallFollow (serPort, bump, dist)
%
% global side;
%
% WALL_CLOSE = 200;
% WALL_FAR = WALL_CLOSE + 200;
% SHARP_TURN = 0.5;
% SLIGHT_TURN = 0.25;
% SPEED = 0.20;
%
% if bump
%     SetFwdVelAngVelCreate(serPort, 0.0, SHARP_TURN*side);
% elseif dist > WALL_FAR || isnan(dist)
%     SetFwdVelAngVelCreate(serPort, SPEED, -side*SLIGHT_TURN);
% elseif dist < WALL_CLOSE
%     SetFwdVelAngVelCreate(serPort, SPEED, side*(dist/WALL_CLOSE*SHARP_TURN));
% else
%     SetFwdVelAngVelCreate(serPort, SPEED, 0.0);
% end
% end
function WallFollow (serPort, startdemo)
if startdemo
    DemoCmdsCreate(serPort, 3); % run the mouse wall follow demo
end
end

% OpenField (behavior state)
%
% Description
%   Moves roomba towards the most open area of the environment as detected
%   by laser readings. If a collision occurs, rotate clockwise.
%
% Inputs
%   serPort - serial port file id for roomba
%   bump - 1 if collision detected
%   left - average distance to the left of the roomba (measured in mm).
%   middle - average distance in front of the roomba (measured in mm).
%   right - average distance to the right of the roomba (measured in mm).
%
% Outputs
%   none
%
function OpenField (serPort, bump, left, middle, right, MINSPEED, MAXSPEED)

global side;

MAXDIST = 3000;
TURN = 0.15;

% find the most open area. speed is proportional to the amount of open space

% first check for a bump
if bump
    SetFwdVelAngVelCreate(serPort, 0.0, 2*TURN*side);
    
    % go straight
elseif (middle > left && middle > right) || isnan(middle)
    if isnan (middle) || middle > MAXDIST
        middle = MAXDIST;
    end
    SetFwdVelAngVelCreate(serPort, max(MINSPEED,(middle/MAXDIST)^2*MAXSPEED), 0.0);
    
    % go left
elseif left > right || isnan(left)
    if isnan (left) || left > MAXDIST
        left = MAXDIST;
    end
    SetFwdVelAngVelCreate(serPort, max(MINSPEED,(left/MAXDIST)^2*MAXSPEED), TURN);
    
    % go right
else
    if isnan (right) || right > MAXDIST
        right = MAXDIST;
    end
    SetFwdVelAngVelCreate(serPort, max(MINSPEED,(right/MAXDIST)^2*MAXSPEED), -1*TURN);
end
end

% ExploreObject (behavior state)
%
% Description
%   Moves roomba towards the area of the environment that changed the most
%
% Inputs
%   serPort - serial port file id for roomba
%   bump - 1 if collision detected
%   left - change to the left of the roomba (measured in mm).
%   middle - change in front of the roomba (measured in mm).
%   right - change to the right of the roomba (measured in mm).
%   mousecnt - number of pixels responding to mouse image;
%   mousectr - centroid of mouse image in horizontal plane
%
% Outputs
%   none
%
function ExploreObject (serPort, bump, objdist, len, left, middle, right)

MINSPEED = 0.10;
TURN = 0.5;


if bump
    SetFwdVelAngVelCreate(serPort, 0.0, TURN);
elseif objdist > 0
    % object on the right
    if objdist < len*3/8
        SetFwdVelAngVelCreate(serPort, MINSPEED, -1*TURN);
        
        % object on the left
    elseif objdist > len*5/8
        SetFwdVelAngVelCreate(serPort, MINSPEED, 1*TURN);
        
        % object in the middle
    else
        SetFwdVelAngVelCreate(serPort, MINSPEED, 0.0);
    end
else
    OpenField (serPort, bump, left, middle, right, 0.10, 0.15);
end

end

% FindHome (behavior state)
%
% Description
%   Set roomba mode to the Dock and Cover demo. Performs random search for
%   the docking station.
%
% Inputs
%   serPort - serial port file id for roomba
%
% Outputs
%   none
%
function FindHome (serPort,startdemo)
if startdemo
    DemoCmdsCreate(serPort, 1); % run the cover and dock demo
end
end

% AtHome (sub-state of FindHome)
%
% Description
%   Set roomba mode to the Dock and Cover demo. Performs random search for
%   the docking station. Set battery level to full by resetting initial
%   battery level.
%
% Inputs
%   serPort - serial port file id for roomba
%
% Outputs
%   none
%
function AtHome (serPort)
global td;
global battery_initial_level;

fwrite(serPort, [143]);
% pause (td);
SetFwdVelAngVelCreate(serPort, 0.0, 0.0);
pause(5);
[batcharge, batcapacity, battery_level] = BatteryChargeReaderRoomba(serPort);
battery_initial_level = battery_level;    % initial battery level
end

% LeaveHome (sub-state of FindHome)
%
% Description
%   Set roomba mode to full. Allows program to regain control of roomba.
%   Back away from docking station, turn ~180 degrees.
%
% Inputs
%   serPort - serial port file id for roomba
%
% Outputs
%   none
%
function LeaveHome (serPort)

fwrite(serPort, [132]);   % go into full mode
SetFwdVelAngVelCreate(serPort, -0.15, 0.0);
pause(2.5);
SetFwdVelAngVelCreate(serPort, 0.0, 0.50);
pause(5);
SetFwdVelAngVelCreate(serPort, 0.0, 0.0);
end

% DockBeam
%
% Description
%   Gets the current state from the roomba's dock beam sensor.
%
% Inputs
%   serPort - serial port file id for roomba
%
% Outputs
%   beam - current status of dock beam sensor. zero means no reading.
%       non-zero denotes reading from nearby docking station.
%
function [beam] = DockBeam(serPort)
% returns the state of the dock beam (green/red buoy, force field, etc)

%Initialize preliminary return values
beam = nan;
try
    
    %Flush Buffer
    N = serPort.BytesAvailable();
    while(N~=0)
        fread(serPort,N);
        N = serPort.BytesAvailable();
    end
    
    warning off
    global td
    
    fwrite(serPort, [142]);  fwrite(serPort,17);
    beam = fread(serPort, 1);
    
    pause(td)
catch
    disp('WARNING:  function did not terminate correctly.  Output may be unreliable.')
end

switch (beam)
    case 248
        % disp ('Red Buoy');
    case 244
        % disp ('Green Buoy');
    case 242
        % disp ('Force Field');
    case 252
        % disp ('Red Buoy and Green Buoy')
    case 250
        % disp ('Red Buoy and Force Field')
    case 246
        % disp ('Green Buoy and Force Field')
    case 254
        % disp ('Red Buoy, Green Buoy and Force Field')
    otherwise
        beam = 0;
        % disp (['No beam: ', num2str(beam)])
end
end

% urg_getdata
%
% Description
%   Read from the URG laser range finder as specified by the URG manual.
%   Plot the ranges in polar coordinates.
%
% Inputs
%   serPort - serial port file id for URG laser range finder
%
% Outputs
%   ranges - reads back ~270 degrees at half a degree resolution. ranges
%   are given in millimeters
%
function ranges = urg_getdata (serPort)

% fprintf(serPort,'MD0044072501000\n'); % scan infinitely
fprintf(serPort,'MD0044072501001\n');   % scan once

[c,n]=fscanf(serPort, '%s',512);
while n < 64
    [c,n]=fscanf(serPort, '%s',512);
    %    [c,n]=fread(serPort,[1,64],'char');
    %     disp(c)
end

inx = 1;
% disp(['n=', num2str(n), ' c=', num2str(size(c,2)), ' inx=' num2str(inx)])

while n > 1
    dat(inx:inx+size(c,2)-2)=c(1:size(c,2)-1);
    inx = inx + size(c,2)-1;
    %     disp(['n=', num2str(n), ' c=', num2str(size(c,2)-1), ' inx=' num2str(inx)])
    [c,n]=fscanf(serPort, '%s',512);
    %    [c,n]=fread(serPort,[1,64],'char');
    % dat(inx:inx+n-3) = c(1:n-2);
    % disp(c)
end

% convert data based on the Hokuyo communication protocol. SCIP 2.0.
inx = 0;
for i=1:size(dat,2)/3
    inx = inx+1;
    ranges(inx)=(dat((i-1)*3+1:i*3)-48)*[2^12 2^6 1]';
end

plotranges = ranges;
plotranges(end+1) = nan;
plotranges(end+1) = 5500;

%     subplot(1,2,1)
%     plot(ranges)

%     subplot(1,2,2)
step_size = 240/682*pi/180;
start = 330/360 * 2*pi;
theta = (1:size(plotranges,2))*step_size + start;
polar(theta,plotranges);
drawnow
end

% urg_getdist
%
% Description
%   Calculate the minimum, average and maximum distances.
%
% Inputs
%   r - array of laser range finder reading
%
% Outputs
%   min - closest range reading
%   avg - average range reading over the array
%   max - farthest range reading
%
function [min, avg, max] = urg_getdist(r)

cnt = 0;
sum = 0;
min = 5600;
max = -1;

for i = 1:size(r,2)
    if r(i) < 5500 && r(i) > 1
        cnt = cnt + 1;
        sum = sum + r(i);
        if r(i) < min
            min = r(i);
        end
        if r(i) > max
            max = r(i);
        end
    end
end
avg = sum/cnt;
end

function y = urg_findobj (r)

x=find(r < 1000);
if mean(x(2:end)-x(1:end-1)) == 1 && size(x,2) > 20 && size(x,2) < 75
    y = mean(x);
else
    y = -1;
end
end

% Roomba Network Initializer
%
%
% Description
%   Instantiates neural network based on principles of the brain's neuromodulatory
%   system. Neurons and connections in the network are specified and
%   initialized.
%
% Inputs
%   none
%
% Outputs
%   none
%
function roomba_net_init (tc_5ht, tc_da, mpfc, ofc)

% state neurons
global STATE_WALL_FOLLOW;
global STATE_OPEN_FIELD;
global STATE_EXPLORE_OBJECT;
global STATE_FIND_HOME;

global EVENT_OBJECT;
global EVENT_LIGHT;
global EVENT_BUMP;
global E;

global N;
global n;
N = 4;
n = zeros (1,N);

% neuromodulators
global NM_DA;
global NM_5HT;
global NM;
global nm;
global tonic;
global TONIC_5HT_DECAY;  % recovery time constant for serotonin
global TONIC_DA_DECAY;  % recovery time constant for dopamine
TONIC_5HT_DECAY = tc_5ht;
TONIC_DA_DECAY = tc_da;

NM_DA =  1;
NM_5HT = 2;

NM = 2;
nm = zeros (1,NM);
tonic = ones (1,NM);
tonic(NM_5HT) = 2;

% ACH/NE neuromodulation
global achne;
achne = zeros (1,E);

% set weights to their initial values
global w_n_n_exc;
global w_n_n_inh;
global w_n_nm;
global w_nm_n;
global w_nm_nm;
global w_achne_n;
global w_e_nm;
global w_e_achne;

% state neuron intrinsic connectivity
w_n_n_exc = 1.0*ones (N,N);
w_n_n_inh = -1.0*ones (N,N);
for i=1:N
    w_n_n_exc (i,i) = 0;
    w_n_n_inh (i,i) = 0;
end

% simulate mPFC projection to raphe to control stressful event
% simulate OFC projection to vta to control rewarding event
w_n_nm = zeros(N,NM);
if mpfc
    w_n_nm (STATE_FIND_HOME, NM_5HT) = -1.0;
    w_n_nm (STATE_WALL_FOLLOW, NM_5HT) = -1.0;
end

if ofc
    w_n_nm (STATE_EXPLORE_OBJECT, NM_DA) = -1.0;
    w_n_nm (STATE_OPEN_FIELD, NM_DA) = -1.0;
end

% neuromodulator to state neuron connectivity
w_nm_n = zeros(NM,N);
w_nm_n (NM_5HT, STATE_FIND_HOME) = 5;
w_nm_n (NM_5HT, STATE_WALL_FOLLOW) = 5;
w_nm_n (NM_DA, STATE_EXPLORE_OBJECT) = 5;
w_nm_n (NM_DA, STATE_OPEN_FIELD) = 5;

% neuromodulatory to neuromodulatory connectivity
%   5HT inhibits DA
w_nm_nm = zeros(NM,NM);
w_nm_nm (NM_5HT, NM_DA) = -1.0;   % change this to change opponency

% event neuron to state neuron connectivity
w_achne_n = 3*ones(E,N);

% event neuron to neuromodulator  connectivity
w_e_nm = zeros(E,NM);
w_e_nm (EVENT_OBJECT, NM_DA) = 0.5;
w_e_nm (EVENT_LIGHT, NM_5HT) = 0.5;
w_e_nm (EVENT_BUMP, NM_5HT) = 0.5;    % risk averse behavior
w_e_nm (EVENT_BUMP, NM_DA) = 0.5;   % risk taking behavior

% event neuron to neuromodulator  connectivity
w_e_achne = ones(1,E);

end

% Roomba Network Cycle
%
% Description
%   Calculates the neural activity and synaptic plasticity for one
%   simulation cycle. The neural network is based on principles of the
%   brain's neuromodulatory system.
%
% Inputs
%   e - a binary array of sensory events.
%
% Outputs
%   choice - an index of a state neuron. zero if no action is selected. non-zero if action selected.
%   achne_out - activity of ACh/NE neurons
%   n_out - activity of state neurons
%   nm_out - activity of neuromodulatory neurons
%
function [choice, Iachne, In, Inm] = roomba_net_cycle (e)

ACTION_SELECTION_THRESHOLD = 0.67;

global E;

% parameters for state neuron activation function
global N;
global n;
N_ACT_GAIN = 2;
N_ACT_PERSIST = 0.25;
N_ACT_BASECURRENT = -1.0;
N_ACT_THR = 0.5;
N_ACT_RESET = 0.9;
nprev = n;

% parameters for neuromodulatory neuron activation function and synaptic
% plasticity
global NM;
global nm;
global tonic;
global NM_DA;
global NM_5HT;
nmprev = nm;
NM_ACT_BASECURRENT = -1.0;
NM_ACT_PERSIST = 0.0;  % persistence of synaptic current
NM_STP_GAIN = 1.0;  % facillitating synapse
NM_STP_DECAY = 25;  % recovery time constant
NM_STP_MAX = 2;     % weight value ceiling
NM_5HT_GAIN = 2;
NM_DA_GAIN = 2;    % gain for sigmoid function

% parameters for ACh/NE neuron activation function and synaptic plasticity
global achne;
achneprev = achne;
ACHNE_ACT_GAIN = 10; % gain for sigmoid function
ACHNE_ACT_BASECURRENT = -0.5;
ACHNE_ACT_PERSIST = 0.5;   % persistence of synaptic current
ACHNE_STP_GAIN = 0.25;   % depressing synapse
ACHNE_STP_DECAY = 25;   % recovery time constant
ACHNE_STP_MAX = 1;      % weight value ceiling

global TONIC_5HT_DECAY;  % recovery time constant for serotonin
global TONIC_DA_DECAY;  % recovery time constant for dopamine
TONIC_STP_GAIN = 1.25;  % facillitating synapse
TONIC_STP_MAX = 2;     % weight value ceiling


global w_n_n_inh;
global w_n_nm;
global w_nm_n;
global w_nm_nm;
global w_achne_n;
global w_e_nm;
global w_e_achne;

% calculate cholinergic/noradrenergic neural activity
for i = 1:E
    I = ACHNE_ACT_BASECURRENT + ACHNE_ACT_PERSIST * achneprev(i) + e(i)*w_e_achne(i);
    achne(i) =  activity (I, ACHNE_ACT_GAIN);
    Iachne(i) = I;
end

% calculate neuromodulatory activity
for i = 1:NM
    %     I = NM_ACT_BASECURRENT + NM_ACT_PERSIST * nmprev(i);
    I = NM_ACT_BASECURRENT + NM_ACT_PERSIST * nmprev(i) + tonic(i);
    for j = 1:E
        I = I + e(j)*w_e_nm(j,i);
    end
    for j = 1:NM
        I = I + nmprev(j)*w_nm_nm(j,i);
    end
    
    for j = 1:N
        I = I + nprev(j)*w_n_nm(j,i);
    end
    
    if i == NM_DA
        nm(i) =  activity (I, NM_DA_GAIN);
    else
        nm(i) =  activity (I, NM_5HT_GAIN);
    end
    Inm(i) = I;
end

% calculate state neural activity
for i = 1:N
    I = N_ACT_BASECURRENT+rand+N_ACT_PERSIST * nprev(i);
    
    % intrinsic synaptic input
    for j = 1:N
        I = I + sum(achneprev) * nprev(j) * w_n_n_inh(j,i);
    end
    
    % event synaptic input
    for j = 1:NM
        I = I + sum(achneprev) * nmprev(j) * w_nm_n(j,i);
    end
    
    %     for j = 1:E
    %         I = I + achneprev(j) * w_achne_n(j,i);
    %     end
    n(i) = activity (I, N_ACT_GAIN);
    In(i) = I;
end

% update plastic weights with short-term plasticity rule. a spike occurs
% when an event occurs.
for i = 1:E
    w_e_achne(i) = stp (w_e_achne(i), ACHNE_STP_GAIN, ACHNE_STP_DECAY, ACHNE_STP_MAX, e(i) > 0.5);
end

% for i = 1:E
%     for j = 1:NM
%         if w_e_nm (i,j) > 0
%             w_e_nm (i,j) = stp (w_e_nm (i,j), NM_STP_GAIN, NM_STP_DECAY, NM_STP_MAX, e(i) > 0.5);
%         end
%     end
% end

for i = 1:E
    if w_e_nm (i,NM_5HT) > 0
        tonic(NM_5HT) = stp (tonic(NM_5HT), TONIC_STP_GAIN, TONIC_5HT_DECAY, TONIC_STP_MAX, achne(i) > 0.5);
    end
    if w_e_nm (i,NM_DA) > 0
        tonic(NM_DA) = stp (tonic(NM_DA), TONIC_STP_GAIN, TONIC_DA_DECAY, TONIC_STP_MAX, achne(i) > 0.5);
    end
end

% find most active state neuron. perform action selection if activity is
% above threshold
[y,i] = max(n);
if y > ACTION_SELECTION_THRESHOLD
    choice = i;
else
    choice = 0;
end

end

% stp - Short-Term Plasticity
%
% Description
%   Simple version of short-term plasticity rule
%
% Inputs
%   xin - current weight value
%   p - amount to increase or decrease weight if there is a spike
%   tau - recovery time constant.
%   max - maximum weight value
%   spk - 1 if spike occurred.
%
% Outputs
%   x - new weight value
%
function x = stp (xin, p, tau, max, spk)

if spk
    x = p*xin;
else
    x = xin + (1-xin)/tau;
end
x = min(max,x);
end

% activity - sigmoid activation function
%
% Description
%   Sigmoid activation function for rate neuron
%
% Inputs
%   I - synaptic input
%   g - gain or slope of sigmoid curve
%
% Outputs
%   s - activity of neuron between 0 and 1
%
function s = activity (I, g)

s = 1/(1+exp(-g*I));
end

% findMouse - find MightyMouse image
%
% Description
%   Finds red background and cape of Mighty Mouse image.
%   Assuming red = activation on R channel and lack of activation on
%   G channel.
%
% Inputs
%   frame - RGB camera frame or image
%
% Outputs
%   pixelcnt - number of pixels that correspond to red.
%   center - centroid of activity along the horizontal plane.
function [pixelcnt, center] = findMouse(frame)

smiley = bwareaopen(frame(:,:,1) > 32 & frame(:,:,2) < 32,10);
pixelcnt = sum(sum(smiley));
center=sum([1:160] .* sum(smiley,1))/sum(sum(smiley));
end
