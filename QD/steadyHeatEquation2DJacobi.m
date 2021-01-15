% steadyHeatEquation2DJacobi.m
% Zijun Fang
% CID: 01811420
% -------------------
% Solve the 2D poisson equation
%
% -(T_xx + T_yy) = q(x,y)
%
% With the boundary conditions on the edge sides, which modelling in the L-
% shape plate with each side grid point in Nx and Ny; to modfy the unit
% length of gird step of finite element mesh to simulate the heat in
% different solving methods.
%
% Solve this question in Jacobi method.

% Clean up
clear
close all

% ----------------------------------------------------------------------- %
% Setup
% ----------------------------------------------------------------------- %

%% Set the dimension of matrix A
Nx = 5;
Ny = 10;
a = 1;
dx = a / (Nx+1);
dy = a / (Ny+1);
A = zeros(Nx*(2*Ny+1)+(Nx+1)*Ny);
c = (dx/dy)^2;
c1 = 2*(1+c);
c2 = c;
c3 = 1;
n = 1:Nx;
m = (Nx*(2*Ny+1)+(Nx+1)*Ny);

%% Matrix A setting
% General diag matrix settings

for i = 1:(Nx*(2*Ny+1)+(Nx+1)*Ny)
    A(i,i) = -c1;
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


%% Setting the vector b which based on function q(x,y)
% Setting the calculation loop of output temp results
b = zeros (m,1);
k = 1;

for i = dx : dx : 2*a-dx
    for j = dy : dy : 2*a-dy
        if i >= a-dx && j <= a+dy
        else
            b(k) = -dx^2 * qxy(i,j,a);
            k = k+1;
        end
    end
end

           

%% Jacobi method simulation and Temperature output
% Form the matrix A = L + U + D
% D : diagonal matrix
% L : lower matrix
% U : upper matrix
D = diag ( diag (A) );
L = tril (A) - D;
U = triu (A) - D;


% Jacobian matrix function : Tx_u = P * b - Q * Tx
% u : time step number
P = D^-1;
Q = (P) * (L+U);

% Temperature output results
Tx = zeros (m,1);
Tx_u = zeros (m,1);

tol_min = 0.0001;
iter_max = 10^4;

for iter = 0 : iter_max
    Tx = P * b - Q * Tx;
    plotT(Nx, Ny, dx, dy, Tx);
    tol = norm( (Tx_u-Tx), 2)
    
    if tol > 100
        disp(' This matrix forming in weak solution.');
        break;
    elseif tol < tol_min
        disp(' This matrix forming in possible solution.');
        break;
    end
    Tx_u = Tx;
    
end


%% Function assign T and plot contour
function plotT(Nx, Ny, dx, dy, Tx)

[X,Y]=meshgrid(0:dy:2,0:dx:2);
T = zeros(2*Nx+3,2*Ny+3);

k=1;
for i = 2 : 2*Nx+2
    for j = 2 : 2*Ny+2
        if i > Nx+1 && j <= Ny+2
        else
        T(i , j) = Tx(k);
        k = k+1;
        end
    end
end

surfc(Y,X,T)
title(['Jacobi Method Heat Equation with:', ' Nx=', num2str(Nx), ', Ny=', num2str(Ny)])
xlabel('Nx-grid');
ylabel('Ny-grid');
zlabel('Temperature');
end
%% Function generating q(x,y)

% q = @(x,y) (100*cos(pi*x) * sin(1.5*pi*y)) * (exp(-0.5*((x-1)^2 +
% (y-1)^2)))

function q = qxy(x,y,a) 
q = (100*cos(pi*x) * sin(1.5*pi*y)) * (exp(-0.5*((x-a)^2 + (y-a)^2)));
end


















