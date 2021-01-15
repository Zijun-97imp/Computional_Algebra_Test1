function A = diffusionMatrix2D(Nx,Ny);
% diffusionMatrix2D.m
% Zijun Fang
% CID: 01811420
% -----------------------------
% solve the 2D diffusion matrix
%
%
% with Dirichlet boundary conditions on the unit square using finite
% differences. The resulting linear system is solved using a direct method.

% Clean up
clear
close all

% ----------------------------------------------------------------------- %
% Setup
% ----------------------------------------------------------------------- %
%% Dimensional settings
Nx = 3;
Ny = 4;
a = 1;
dx = a / Nx;
dy = a / Ny;
A = zeros(Nx*(2*Ny+1)+(Nx+1)*Ny);
c = (dx / dy)^2;
c1 = 2*(1+c);
c2 = c;
c3 = 1;
n = 1:Nx;


%% Matrix Settings
% General diag matrix settings

for i = 1:(Nx*(2*Ny+1)+(Nx+1)*Ny)
    A(i,i) = c1;
end
for i = 2:(Nx*(2*Ny+1)+(Nx+1)*Ny)
    A(i,i-1) = c2;
    A(i-1,i) = c2;
end
for i = (2*Ny+2):(Nx*(2*Ny+1))
    A(i-(2*Ny+1),i) = c3;
    A(i,i-(2*Ny+1)) = c3;
end

for i = (Nx*(2*Ny+1)+1+Ny):(Nx*(2*Ny+1)+(Nx+1)*Ny)
    A(i-Ny,i) = c3;
    A(i,i-Ny) = c3;
end

% Specific matrix settings
for i = n*(2*Ny+1)
    A(i,i+1) = 0;
    A(i+1,i) = 0;
end

for i = n*Ny+(Nx*(2*Ny+1))
    A(i,i+1) = 0;
    A(i+1,i) = 0;
end

A = sparse(A);

spy(A)
        














