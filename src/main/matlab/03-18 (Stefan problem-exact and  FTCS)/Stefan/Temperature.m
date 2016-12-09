%Conversion from conservative variable q (internal energy) to the primitive
%variable T (temperature)
function T=Temperature(q)
global kL kR hc cL cR Tc TR TL rhoL rhoR lambdaL lambdaR
if (q<0)
    %pure ice
    T = q/(rhoL*cL);
elseif (q>=0 && q<=rhoR*hc)
    %phase transition region
    T=0;
else
    %pure liquid water
    T=(q-rhoR*hc)/(rhoR*cR);
end

