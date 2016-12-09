%Thomas agorithm for the solution of tridiagonal system
%Input:
%   a=vector of lower diagonal elements
%   b= vector of diagonal elements
%   c=vector of upper diagonal elements
%Output
% x= solution vector

%Function name:Thomas
%File name:Thomas.m (compulsory)

function x=Thomas(a,b,c,d)  %we need to know how many incognite, posso aggoungere un input oppure
N= length(d); %get the number of elemnts (problem size)
%Part 1 :forward elimination
c(1)=c(1)/b(1);
d(1)=d(1)/b(1);
for i=2:N
    gamma=1/(b(i)-c(i-1)*a(i));
    c(i)=c(i)*gamma;
    d(i)=(d(i)-a(i)*d(i-1))*gamma;
    
end
%part 2: back substitution
x(N)=d(N);
for i=N-1:-1:1      %up with negative step
    x(i)=d(i)-c(i)*x(i+1);
end
x=x(:); %force matlab to produce x as a column vector,