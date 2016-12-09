% internal energy as a function of temperature 
% using a piecewise linear regularization of the
% jump due to the latent heat 
function y=q(T)
global kL kR hc cL cR rhoL rhoR lambdaL lambdaR Tc TL TR epsilon 
if(T<=Tc-epsilon) 
    % pure ice 
    y = rhoL*cL*(T-Tc); 
elseif(T>=Tc+epsilon)
    % pure liquid water 
    y = rhoR*cR*(T-Tc) + rhoR*hc; 
else
    % phase transition region (due to the regularization) 
    dQdT=(rhoR*cR*(Tc+epsilon-Tc) + rhoR*hc - ...
          rhoL*cL*(Tc-epsilon-Tc) )/(2*epsilon); 
    y = rhoL*cL*(Tc-epsilon-Tc) + dQdT*(T-(Tc-epsilon));   
end
