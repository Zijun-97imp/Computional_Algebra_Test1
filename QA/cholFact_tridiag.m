function R = cholFact_tridiag(A)

% R = cholFact(A) returns the upper triangular Cholesky factor of a
%     positive semidefinite symmetric matrix A. This function only operates
%     on the upper-triangular part of A and implements the CHolesky
%     algorithm.
%
% Zijun Fang
% CID: 01811420

% Check inputs
[m, n] = size(A);
assert(m==n, 'Incorrect input size: A must be n-by-n')

% Zero all entries below the main diagonal
R = triu(A);

% Minimum location to calculate the loop in j-direction
% pmj = min(j+2,n);
% Minimum location to calculate the loop in i-direction
% pmi = min(i+1,n);
% Operate on the upper triangular part using Cholesky algorithm (see
% Algorithm 4.11 of the lecture notes)
for i = 1:n
   if R(i,i) > 0
      for j = (i+1):n
        pmj = min(j+2,n);
        R(j,j:min(j+2,n)) = R(j,j:pmj) - R(i,j:pmj).*( R(i,j)/R(i,i) );
      end
      pmi = min(i+1,n);
      R(i,i:pmi) = R(i,i:pmi)./sqrt(R(i,i));
   else
       error('A is not positive definite!')
   end
end