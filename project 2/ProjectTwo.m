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
Ts = 1/fs;
ep = 0.1;
LbAtck = ep ^ (Ts/(2/1000));
LbRels = ep ^ (Ts/(1/100));

%% Part C

fo = [0.3 0.6 1]*1000;
Ao = [2 4 0.5];
to = [0 25 50]/1000;
Smp = 0.025*fs;     % number of samples per 25ms section
time = @(n,t0)  n/(Smp-1) *0.025 + t0;
x = zeros(1,Smp*3);

for i =1:3
    for n = 1:Smp
        t = time(n-1,to(i));
        x(n + 200*(i-1)) = X(t,Ao(i),fo(i));
    end
end


% Mean absolute values
x1MA = AbsAvg(Ao(1));
x2MA = AbsAvg(Ao(2));
x3MA = AbsAvg(Ao(3));

n = 1:length(x);
figure;
plot(n/8,x(n));

%% Part D

p = 1/3;
c0 = 1;
D = 0;
L = 50;             %NEED TO FIND PROPERLY 

MvAvgBuf = zeros(1,L);
pt = @(i) mod(i,L-1) + 1;
MvAvgSum = 0;
MvAvgOld = 0;
cPrev = 0;
c = zeros(1,length(x));
g = zeros(1,length(x));
G = zeros(1,length(x));
y = zeros(1,length(x));

% Level detector

cN = @(xn,cPrev) (LbAtck* cPrev + (1-LbAtck)*abs(xn)).*(abs(xn) >= cPrev) ...
    + (LbRels*cPrev + (1 - LbRels).*abs(xn)).*(abs(xn) < cPrev);

% Gain Processor

gCom = @(c) (c/c0)^(p-1)*(c>= c0) + 1*(c <=c0);

gExp = @(c) (c/c0)^(p-1)*(c <=c0) + 1*(c>= c0);

% Since we are compressing the signal I will use gCom
for i = 1:length(x)

    xn = x(i);
    c(i) = cN(xn,cPrev);
    cPrev = c(i);
    g(i) = gCom(c(i));
    MvAvgOld = MvAvgBuf(pt(i-1));
    MvAvgBuf(pt(i-1)) = g(i);
    MvAvgSum = MvAvgSum + g(i) - MvAvgOld;
    G(i) = MvAvgSum/L;
    y(i) = G(i)*xn;
    
end

figure;
plot(n,y)

figure;
plot(n,c)

figure;
plot(n,g)

figure;
plot(n,G)

%% Part E Limiter

p = 1/15;
c0 = 2;
D = 0;
L = 40;             %NEED TO FIND PROPERLY 

MvAvgBuf = zeros(1,L);
pt = @(i) mod(i,L-1) + 1;
MvAvgSum = 0;
MvAvgOld = 0;
cPrev = 0;
c = zeros(1,length(x));
g = zeros(1,length(x));

% Level detector

cN = @(xn,cPrev) (LbAtck* cPrev + (1-LbAtck)*abs(xn)).*(abs(xn) >= cPrev) ...
    + (LbRels*cPrev + (1 - LbRels).*abs(xn)).*(abs(xn) < cPrev);

% Gain Processor

gCom = @(c) (c/c0)^(p-1)*(c>= c0) + 1*(c <=c0);

gExp = @(c) (c/c0)^(p-1)*(c <=c0) + 1*(c>= c0);

% Since we are compressing the signal I will use gCom
for i = 1:length(x)

    xn = x(i);
    c(i) = cN(xn,cPrev);
    cPrev = c(i);
    g(i) = gCom(c(i));
    MvAvgOld = MvAvgBuf(pt(i-1));
    MvAvgBuf(pt(i-1)) = g(i);
    MvAvgSum = MvAvgSum + g(i) - MvAvgOld;
    G(i) = MvAvgSum/L;
    y(i) = G(i)*xn;
    
end

figure;
plot(n,y)

figure;
plot(n,c)

figure;
plot(n,g)

figure;
plot(n,G)

%% Part F Expander

p = 4;
c0 = 0.95;
D = 0;
L = 40;             %NEED TO FIND PROPERLY 

MvAvgBuf = zeros(1,L);
pt = @(i) mod(i,L-1) + 1;
MvAvgSum = 0;
MvAvgOld = 0;
cPrev = 0;
c = zeros(1,length(x));
g = zeros(1,length(x));

% Level detector

cN = @(xn,cPrev) (LbAtck* cPrev + (1-LbAtck)*abs(xn)).*(abs(xn) >= cPrev) ...
    + (LbRels*cPrev + (1 - LbRels).*abs(xn)).*(abs(xn) < cPrev);

% Gain Processor

gCom = @(c) (c/c0)^(p-1)*(c>= c0) + 1*(c <=c0);

gExp = @(c) (c/c0)^(p-1)*(c <=c0) + 1*(c>= c0);

% Since we are compressing the signal I will use gCom
for i = 1:length(x)

    xn = x(i);
    c(i) = cN(xn,cPrev);
    cPrev = c(i);
    g(i) = gExp(c(i));
    MvAvgOld = MvAvgBuf(pt(i-1));
    MvAvgBuf(pt(i-1)) = g(i);
    MvAvgSum = MvAvgSum + g(i) - MvAvgOld;
    G(i) = MvAvgSum/L;
    y(i) = G(i)*xn;
    
end

figure;
plot(n*3/25,y)

figure;
plot(n*3/25,c)

figure;
plot(n*3/25,g)

figure;
plot(n*3/25,G)



%% Part G Noise Gate: f1 & f3

p = 12;
c0 = 1.9;
D = 0;
L = 40;             %NEED TO FIND PROPERLY 

MvAvgBuf = zeros(1,L);
pt = @(i) mod(i,L-1) + 1;
MvAvgSum = 0;
MvAvgOld = 0;
cPrev = 0;
c = zeros(1,length(x));
g = zeros(1,length(x));

% Level detector

cN = @(xn,cPrev) (LbAtck* cPrev + (1-LbAtck)*abs(xn)).*(abs(xn) >= cPrev) ...
    + (LbRels*cPrev + (1 - LbRels).*abs(xn)).*(abs(xn) < cPrev);

% Gain Processor

gCom = @(c) (c/c0)^(p-1)*(c>= c0) + 1*(c <=c0);

gExp = @(c) (c/c0)^(p-1)*(c <=c0) + 1*(c>= c0);

% Since we are compressing the signal I will use gCom
for i = 1:length(x)

    xn = x(i);
    c(i) = cN(xn,cPrev);
    cPrev = c(i);
    g(i) = gExp(c(i));
    MvAvgOld = MvAvgBuf(pt(i-1));
    MvAvgBuf(pt(i-1)) = g(i);
    MvAvgSum = MvAvgSum + g(i) - MvAvgOld;
    G(i) = MvAvgSum/L;
    y(i) = G(i)*xn;
    
end

figure;
plot(n,y)

figure;
plot(n,c)

figure;
plot(n,g)

figure;
plot(n,G)



%% Part G Noise Gate:  f3

p = 11;
c0 = 1.25;
D = 0;
L = 40;             %NEED TO FIND PROPERLY 

MvAvgBuf = zeros(1,L);
pt = @(i) mod(i,L-1) + 1;
MvAvgSum = 0;
MvAvgOld = 0;
cPrev = 0;
c = zeros(1,length(x));
g = zeros(1,length(x));

% Level detector

cN = @(xn,cPrev) (LbAtck* cPrev + (1-LbAtck)*abs(xn)).*(abs(xn) >= cPrev) ...
    + (LbRels*cPrev + (1 - LbRels).*abs(xn)).*(abs(xn) < cPrev);

% Gain Processor

gCom = @(c) (c/c0)^(p-1)*(c>= c0) + 1*(c <=c0);

gExp = @(c) (c/c0)^(p-1)*(c <=c0) + 1*(c>= c0);

% Since we are compressing the signal I will use gCom
for i = 1:length(x)

    xn = x(i);
    c(i) = cN(xn,cPrev);
    cPrev = c(i);
    g(i) = gExp(c(i));
    MvAvgOld = MvAvgBuf(pt(i-1));
    MvAvgBuf(pt(i-1)) = g(i);
    MvAvgSum = MvAvgSum + g(i) - MvAvgOld;
    G(i) = MvAvgSum/L;
    y(i) = G(i)*xn;
    
end

figure;
plot(n,y)

figure;
plot(n,c)

figure;
plot(n,g)

figure;
plot(n,G)



    
