%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Jeff Krichmar
%  Psych 268A - Fall 2015
%  Homework 2
%  Problem 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Em = zeros(1,100);

% concentrations and other parameters taken from lecture slides
c= 58;
pNaK = 0.02;
pClK = 1;
Ki = 125;
Ko = 5;
Nai = 12;
Nao = 120;
Cli = 5;
Clo = 125;

T = 100;

for i=1:100
    if i > 44 && i < 56
        pNaK = 20;
    else 
        pNaK = 0.02;
    end
    
    % note that matlab log is the natural log. use log10 instead.
    Em(i) = c*log10((Ko + pNaK*Nao + pClK*Cli)/(Ki + pNaK*Nai + pClK*Clo));
end

plot(1:T,Em);
title('Relation Between Na+ Permeability and Membrane Potential');
xlabel('Time');
ylabel('Membrane Potential');