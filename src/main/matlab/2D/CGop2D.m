%matrix-free conjugate method for the solution of A*x=b with A symmetric positive definite
%Input: (we don't need the matrix as input)
%       b = right hand side
%output: 
%       x = solution vector
%the user has to provide a function "matop2D(x)" wich computes the matrix-vector product A*x

function x=CGop2D(b)
x = b;                %initial guess MATRIX!
N = numel(b);         %get the problem size 
r = b-matop2D(x);     %residual MATRIX!
p = r;                %initial search direction MATRIX!
alpha = sum(sum(r.*r));  %square of r (scalar product) %or r'*r; 
for k=1:4*N           %2* or 4*N to have more precision
    if (sqrt(alpha)<1e-14)   %nominal machine operation is 1e-16 so we use something smaller
        %norm of the residual is sufficently small
        break
    end
    v = matop2D(p);   % MATRIX!
    lambda = alpha/sum(sum(p.*v)); % MATRIX!
    x = x + lambda*p; % MATRIX!
    r = r - lambda*v; % MATRIX!
    alphak = alpha;
    alpha = sum(sum(r.*r));
    p = r + alpha/alphak*p; % MATRIX!
end


