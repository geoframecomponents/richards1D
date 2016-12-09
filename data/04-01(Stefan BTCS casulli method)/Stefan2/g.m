%function g(gamma) for the Stefan problem
function res = g(gamma) 
global kL kR hc cL cR Tc TR TL rhoL rhoR lambdaL lambdaR
%
res = gamma*sqrt(kL)*hc*rhoR- ...
      lambdaL*(Tc-TL)/erf( gamma            )*exp(-gamma^2      )/sqrt(pi*kL)- ...
      lambdaR*(Tc-TR)/erfc(gamma*sqrt(kL/kR))*exp(-gamma^2*kL/kR)/sqrt(pi*kR);
