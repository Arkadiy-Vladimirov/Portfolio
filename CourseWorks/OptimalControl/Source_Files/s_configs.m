clear;clc;
A = @A1t; B = @B1t; f = @f1t;
Pconf.a = 3; Pconf.b = 3; Pconf.c = 5;
Pconf.alpha = 1; Pconf.beta = 2;
r = 3.5; pnts = [12 10; 10 15; 15 9; 17 15];

    params.N = 100;     %size of initial conditions sample;
    params.dt = 0.001;   %time accuracy
    params.maxT = 0.7;      %max time 
    params.phi = [0 2*pi];  %interval of conjugate directions sample
save configs/config1.mat

clear;clc;
A = @A1t; B = @B1t; f = @f1t;
Pconf.a = 3; Pconf.b = 3; Pconf.c = 5;
Pconf.alpha = 1; Pconf.beta = 2;
r = 3.5; pnts = [12 12; 10 17; 15 11; 17 17];

    params.N = 100;     %size of initial conditions sample;
    params.dt = 0.001;   %time accuracy
    params.maxT = 1.4;      %max time 
    params.phi = [0 2*pi];  %interval of conjugate directions sample
save configs/config2.mat

clear;clc;
A = @A3t; B = @B3t; f = @f3t;
Pconf.a = 0.5; Pconf.b = 3; Pconf.c = 2;
Pconf.alpha = 1; Pconf.beta = 2;
r = 3; pnts = [10 10; 10 15; 15 10; 15 15];

    params.N = 75;
    params.dt = 0.001;
    params.maxT = 0.4;
    params.phi = [0 2*pi];
save configs/config3.mat

%% functions conf1
function A = A1t(t)
    A(1,1) = sin(t); A(1,2) = 2;
    A(2,1) = exp(-t); A(2,2) = cos(t);
end

function B = B1t(t)
    B(1,1) = 1; B(1,2) = exp(2*t);
    B(2,1) = t.^3+2*t-1; B(2,2) = t;
end

function f = f1t(t)
    f = zeros(2,1);
    f(1) = 2*t;
    f(2) = -t;
end

%% functions conf3
function A = A3t(t)
    A(1,1) = 2; A(1,2) = 2;
    A(2,1) = 3; A(2,2) = 0;
end

function B = B3t(t)
    B(1,1) = 1; B(1,2) = 6;
    B(2,1) = 20; B(2,2) = 3;
end

function f = f3t(t)
    f = zeros(2,1);
    f(1) = 8;
    f(2) = 7;
end