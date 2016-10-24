
function actor_critic_maze

clear all
close all

epa = 0.9;
epc = 0.9;
beta = 0.5;

act(1:100,1:2) = 0;
mtbl (1:100,1:2,1:3) = 0;
wgt(1:100,1:3) = 0;

for i = 1:10
    rand('seed',i);
    [a, m, w] = actorcritic_maze(epa, epc, beta, true);
    act = act + a;
    mtbl = mtbl + m;
    wgt = wgt + w;
end

subplot (1,3,1)
plot (wgt(:,1)/10);
axis ([1 100 0 5])
axis square
title ('w(A)')
xlabel('trial')
ylabel ('w')

subplot (1,3,2)
plot (wgt(:,2)/10);
axis ([1 100 0 5])
axis square
title ('w(B)')
xlabel('trial')
ylabel ('w')

subplot (1,3,3)
plot (wgt(:,3)/10);
axis ([1 100 0 5])
axis square
title ('w(C)')
xlabel('trial')
ylabel ('w')

act(1:100,1:2) = 0;
mtbl (1:100,1:2,1:3) = 0;
wgt(1:100,1:3) = 0;

for i = 1:10
    rand('seed',i);
    [a, m, w] = actorcritic_maze(epa, epc, beta, false);
    act = act + a;
    mtbl = mtbl + m;
    wgt = wgt + w;
end

mtbl = mtbl / 10;
for i = 1:100
    pA(i,:) = mysoftmax(mtbl(i,:,1)',beta);
    pB(i,:) = mysoftmax(mtbl(i,:,2)',beta);
    pC(i,:) = mysoftmax(mtbl(i,:,3)',beta);
end

figure;
subplot (1,3,1)
plot (pA(:,1));
axis ([1 100 0 1])
axis square
title ('u=A')
xlabel('trial')
ylabel ('p[L;u]')

subplot (1,3,2)
plot (pB(:,1));
axis ([1 100 0 1])
axis square
title ('u=B')
xlabel('trial')
ylabel ('p[L;u]')

subplot (1,3,3)
plot (pC(:,1));
axis ([1 100 0 1])
axis square
title ('u=C')
xlabel('trial')
ylabel ('p[L;u]')

end

function [atrial, mtrial, wtrial] = actorcritic_maze (beta, epsilon_actor, epsilon_critic, policy_evaluation)
A = 1;
B = 2;
C = 3;
D = 4;
E = 5;
F = 6;
G = 7;

L = 1;
R = 2;

M(1:2,1:3) = 0;
w(1:3) = 0;
v(1:7) = 0;

pos(1) = A; % start at A
v(A) = w(A);

for trial = 1:100
    for t = 1:2

        p = mysoftmax (M(:,pos(t)), beta);
        a(t) = (random('unif', 0, 1) > p(1)) + 1;

        if pos(t) == A && a(t) == L
            pos (t+1) = B;
            r = 0;
        elseif pos(t) == A && a(t) == R
            pos (t+1) = C;
            r = 0;
        elseif pos(t) == B && a(t) == L
            pos (t+1) = D;
            r = 0;
        elseif pos(t) == B && a(t) == R
            pos (t+1) = E;
            r = 5;
        elseif pos(t) == C && a(t) == L
            pos (t+1) = F;
            r = 2;
        else % pos(t) == C && a(t) == R
            pos (t+1) = G;
            r = 0;
        end

        delta = critic (r, v(pos(t+1)), v(pos(t)));
        w (pos(t)) = w(pos(t)) + epsilon_critic * delta;
        v(pos(t)) = w(pos(t));

        if ~policy_evaluation
            M(:,pos(t)) = M(:,pos(t)) + epsilon_actor .* actor(a(t), p) .* delta;
        end

        atrial(trial,t) = a(t);
    end
    wtrial(trial, :) = w;
    mtrial(trial,:,:) = M;
end
end

function y = actor (a, p)

if a == 1
    y(1) = 1 - p(1);
    y(2) = -p(2);
else
    y(1) = -p(1);
    y(2) = 1 - p(2);
end

y = y';

end

function y = critic (r, vnext, v)
y = r + vnext - v;
end

function y = mysoftmax (m, beta)

for i = 1:size(m,1)
    y(i) = exp(beta*m(i))/sum(exp(beta*m));
end

end

