function filterArray = MakeFilterArray

size = 15;
sigmaDir = 3;
sigmaSpeed = 1.16;
prefDir = linspace(0,2*pi-pi/4,8);
prefSpeed = [2,4,8,16,32];
nFilters = length(prefDir)*length(prefSpeed);


[X,Y] = meshgrid(1:size);
filterArray = cell(nFilters,1);
for d = 1:length(prefDir)
    d0 = prefDir(d);
    
    for s = 1:length(prefSpeed)
        s0 = prefSpeed(s);
        
        dMt = zeros(size); 
        sMt = zeros(size);
        for x = 1:size
            
            for y = 1:size
                
                dir = atan2(Y(x,y),X(x,y));
                speed = sqrt(X(x,y).^2+Y(x,y).^2);
                
                dMt(x,y) = exp(sigmaDir*cos(dir-d0)-1);
                sMt(x,y) = exp(-log((speed+s0)/(s0+.33)).^2./(2*sigmaSpeed));
            end
            
        end
        filterArray{s,d} = dMt.*sMt;
        figure;
        mesh(filterArray{s,d});
    end
    
    
end


end