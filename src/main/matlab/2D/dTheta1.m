function y=dTheta1(psi)
global alpha thetas thetar n m Ks psic
if (psi<=psic)
    %left of critical value, take the original derivative
    y = dTheta(psi);
else
    %on the right of the critical value, keep the maximum derivative
    y = dTheta(psic);
end
