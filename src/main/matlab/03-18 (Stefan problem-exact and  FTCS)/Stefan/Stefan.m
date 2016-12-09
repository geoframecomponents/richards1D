clear all 
close all 
clc
global kL kR hc cL cR Tc TR TL rhoL rhoR lambdaL lambdaR
%Phisical parameters
lambdaL = 2.09;          %heat conductivity of ice
lambdaR = 0.6;           %heat conductivity water
hc      = 334e3;         %latent heat
rhoL    = 917;           %density of ice
rhoR    = 1000;          %density of water
cL      = 2108;          %heat capacity of ice
cR      = 4187;          %heat capacity of water
kL = lambdaL/(rhoL*cL);
kR = lambdaR/(rhoR*cR);
Tc      = 0;             %critical melting temperature
TL      = -10;           %air temperature
TR      = 5;             %initial lake temperature
%Domain
xL = 0;
xR = 2;                  %bottom of lake
IMAX = 100;              %number of control volumes
dx = (xR-xL)/IMAX;
x = linspace(xL+dx/2,xR-dx/2,IMAX);
day = 24*3600;           %no. of seconds of day
%plot the exact solution after a certain time
time=2*day;
gamma = Newton(1)     %the argument of Newton is the initial guess

G1 = (Tc-TL)/erf(gamma);
G2 = (Tc-TR)/erfc(gamma*sqrt(kL/kR));
s = gamma*2*sqrt(kL*time);     %interface position
xe = linspace(xL,xR,10*IMAX);  %fine grid to plot exact solution
for i=1:10*IMAX
    if(xe(i)<s)
        %left branch of the solution
        Te(i) = TL + G1*erf(xe(i)/2/sqrt(kL*time));
    else
        %right branch of the solution
        Te(i) = TR + G2*erfc(xe(i)/2/sqrt(kR*time));
    end
end
plot(xe,Te,'r-')    %plot exact solution
hold on
plot(s,Tc,'r+')     %plot the interface (s and critical temperature)

%Afeter two days of -10° we get the interface in "+" (we can plot s:
%s=0.1375--13cm). We have used the true parameters for water,we have not
%jet taking account the heat transfer from air to water, we assume that the hit
%transfer is perfect so we are really able to have exactly -10°.






