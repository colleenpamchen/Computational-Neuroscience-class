%% Midterm Assignment

% Make gradient 
% freq = 8;
% angle = 90; 
% Period= 128;
% xx = freq * cosd(angle);
% yy = freq * sind(angle);
% [X Y] = meshgrid(1:Period);
% grating = cos(2*pi*((xx*X/Period) +(yy*Y/Period)));
% imshow(grating); 
% 
% C = 50; % contrast 
% 0.5* C*X

% Make Stimulus:
prefDir = [-45:45:135]; % preferred orientation (degrees) 
contrast = [0;12;25;50;75]; % ??????????? HOW to use contrast with orientation?? 
degree = [ cos(prefDir);sin(prefDir) ];


% how do I create input neurons using cosine turning curves? 
s = -90:1:180; % orientation angles in degrees
smax= prefDir; % preferred orientation 
fr=zeros(length(prefDir),length(s)); 
rmax= 52.14;
sig= 14.73;

for k=1:length(prefDir)
    
y = rmax * exp( (-1/2)* ( ( s- prefDir(k)) / sig ).^2 );
fr(k,:) = y;

% figure
% plot(s,y)
% title('Gaussian tuning curve of avg. firing rate','FontSize',12,'FontWeight','bold')
% xlabel('s: orientation angle in degrees','FontSize',12,'FontWeight','bold') 
% ylabel('f(Hz)','FontSize',12,'FontWeight','bold') 
% legend( num2str(smax(k)) ,'Location','SouthEast');

end
w_excitatory = [];
w_inhibitory = -1*ones(1,length(s));
w =ones(length(prefDir),length(s));

INH = dot(w,fr);

for k=1:length(smax)
    
    inh_fr = INH / ()
    
end


% dot(INH,w_inhibitory);






