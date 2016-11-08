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
ylabel('Sodium Permeability')
title('Permeability and Potential')