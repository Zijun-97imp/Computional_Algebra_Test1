function x = Rsolve_tridiag(R,b)

% x = Rsolve(R,b) solves the upper triangular system R*x=b using backward
% substitution. Only the upper triangular part of R is used.
%
% Zijun Fang
% CID: 01811420

% Preliminary stuff: size of R and preallocate solution (column vector)
% for speed
n = size(R,1);
x = zeros(n,1);

% minimum element point of generate
% posm = min(i+2,n);

% Compute: backward substitution with no error checks
% Implement the backward substitution formula using vector operations.
x(n) =  b(n)/R(n,n);
for i = n-1:-1:1
    posm = min(i+2,n);
    x(i) =  ( b(i) - dot( R(i,i+1:posm), x(i+1:posm) ) ) / R(i,i);
end