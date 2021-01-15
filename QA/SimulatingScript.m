% Question1(c)Script
% function in diffusionatrix1D.m
% Setting general dimensions of simulation

% -----------------------------------------------------------
% Setup
% Matrix and plotting running
% -----------------------------------------------------------
% Zijun Fang
% CID: 01811420

% Clean up
clear
close all


%% General Dimension Settings
% Basic dimension settings
L = 1;
N = 100;
dt = 10^-2;
dx = 1 / (N+1);
sig = dt / (dx^2);

% Calculation and Centre points situation
x = 0 : dx : L;
xcent = 0+dx/2 : dx : L-dx/2;
kcent = (10^-2) + exp(-10 * (xcent));


%% Boundary Condition and Matrix
% Boundary Conditions and Matrices Settings
% diffusionMatrix1D.m for matrix B = (N,kcent)
% for the time step relation: P*u^(n+1) = Q*u^(n)
B = diffusionMatrix1D (N , kcent);
I = eye(N);
u = 1 - cos(4 * pi * x);
P = I + (sig/2)*B;
Q = I - (sig/2)*B;
u = u';

%% Contour Plot and Time step forward situation
% Calculating the (x,u) plot through time pass
% Along the LU Cholesky-algorithm: Q*u^(n+1) = b;
% P*u^n = (RR')*u^n = b;

for T = 0 : dt : 1.5
    switch T
        case 0
            plot(x,u)
            drawnow
            hold on
        case 0.1
            plot(x,u)
            drawnow
            hold on
        case 0.4
            plot(x,u)
            drawnow
            hold on
        case 1.5
            plot(x,u)
            drawnow
            title('1D Diffusion plot (x,u) @t=0, 0.1, 0.4, 1.5')
            xlabel('x-length of the bar')
            ylabel('u-concentration magnitude')
            hold off
    end
    b = Q * u(2:N+1);
    R = cholFact_tridiag(P);
    u(2:N+1) = Rtsolve_tridiag(R,b);
    u(2:N+1) = Rsolve_tridiag(R,u(2:N+1));
end











