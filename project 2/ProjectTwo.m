% Kevin Quizhpi
% DSP Design
% Project 2
% 9/26/17

%% Part A

X = @(t,A,f) A*cos(2*pi*f*t);

MeanSqAvg = @(A) A^2 /2;
AbsAvg = @(A) 2*A /pi;

%% Part B

fs = 8000;
ep = 0.1;
LbAtck = ep ^ 1/(fs*2/1000);
LbRels = ep ^ 1/(fs*1/100);

%% Part C

fo = [0.3 0.6 1]*1000;
Ao = [2 4 0.5];
to = [0 25 50];

for i =1:3
    for n = to(i):((25*i/1000)*fs-1)
        t = n/(25*i/1000);
        x(n+1) = X(t,Ao(i),fo(i));
    end
end

n = 1:length(x);
figure;
plot(n/8,x(n));
    
