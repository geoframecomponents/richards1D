%derivative of the function g(gamma)for the Stefan problem using a finite difference 
%approximation 2order
function res = dg(gamma)
eps = 1e-7;
res = ( g(gamma+eps)- g(gamma-eps) )/(2*eps);