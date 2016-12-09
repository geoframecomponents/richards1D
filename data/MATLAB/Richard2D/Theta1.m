function y = Theta1(psi)
global alpha thetas thetar n m Ks psic
if (psi<=psic)
    y = Thetaf(psi);
else
    y = Thetaf(psic) + dTheta(psic)*(psi-psic);    
end
