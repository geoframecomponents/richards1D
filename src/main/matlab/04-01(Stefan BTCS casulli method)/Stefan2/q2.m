% internal energy as a function of temperature 
% using a piecewise linear regularization of the
% jump due to the latent heat 
function y=q2(T)
global kL kR hc cL cR rhoL rhoR lambdaL lambdaR Tc TL TR epsilon 
y = q1(T)-q(T); 
