clear all 
close all 
clc
global kL kR hc cL cR Tc TR TL rhoL rhoR lambdaL lambdaR
%Phisical parameters
lambdaL = 2.09;         %heat conductivity of ice
lambdaR = 0.6;          %heat conductivity water
hc      = 334e3;        %latent heat
rhoL    = 917;          %density of ice
rhoR    = 1000;         %density of water
cL      = 2108;         %heat capacity of ice
cR      = 4187;         %heat capacity of water
kL      = lambdaL/(rhoL*cL);
kR      = lambdaR/(rhoR*cR);
Tc      = 0;             %critical melting temperature
TL      = -10;           %air temperature
TR      = 5;             %initial lake
%Domain
xL = 0;
xR = 2;                 %bottom of lake
IMAX = 100;             %number of control volumes
dx = (xR-xL)/IMAX;
x = linspace(xL+dx/2,xR-dx/2,IMAX);
day = 24*3600;          %no. of seconds of day
tend = 2*day;           %set the final simulation time
time = 0;               %initial time
gamma = Newton(1)       %the argument of Newton is the initial guess
%set the initial condition (the right temperature-than it start freezing from the left)
for i=1:IMAX
    T(i) = TR;      % >0  (is useless)
    q(i) = hc*rhoR + (T(i)-Tc)*cR*rhoR;  %we need to initialize the energy not the temperature
end
NMAX=100000;
for i=1:NMAX
    dt = 0.45*dx^2/max(kL,kR); %time step restriction (FTCS scheme)
    if(time+dt>tend)
        dt=tend-time;
    end
    
    %compute temperature and lambda (heat conduction coefficient)
    for i=1:IMAX
        T(i)=Temperature(q(i));
        if (T(i)<Tc)
            lambda(i) = lambdaL;
        else
            lambda(i) = lambdaR;
        end
    end
    
    plot(x,T,'o')       %%jumping
    %plot(x,q,'o')      %%energy... 
    title(sprintf('Current time = %f',time))
    drawnow
    if(time>=tend)
        break
    end
    
    %FTCS finte volume scheme
    for i=1:IMAX
        if(i==1)
            lambdap = 0.5*(lambda(i)+lambda(i+1));
            lambdam = lambda(i);
            fp = -lambdap*(T(i+1)-T(i))/dx;
            fm = -lambdam*(T(i)  - TL )/(dx/2);
        elseif(i==IMAX)
            lambdap = lambda(i);
            lambdam = 0.5*(lambda(i)+lambda(i-1));
            fp = -lambdap*(TR  - T(i)  )/(dx/2);
            fm = -lambdam*(T(i)- T(i-1))/dx;
        else
            lambdap = 0.5*(lambda(i)+lambda(i+1));
            lambdam = 0.5*(lambda(i)+lambda(i-1));
            fp = -lambdap*(T(i+1)- T(i)  )/dx;
            fm = -lambdam*(T(i)  - T(i-1))/dx;
        end
        qnew(i) = q(i)-dt/dx*(fp-fm); %finite volume scheme
    end
    time = time+dt;     %advance time
    q = qnew;           %overwrite the old solution             
end


G1=(Tc-TL)/erf(gamma);
G2=(Tc-TR)/erfc(gamma*sqrt(kL/kR));
s = gamma*2*sqrt(kL*time);    %interface position
xe= linspace(xL,xR,10*IMAX);  %fine grid to plot exact sol
for i=1:10*IMAX
    if(xe(i)<s)
        %left branch of the solution
        Te(i)= TL + G1*erf(xe(i)/2/sqrt(kL*time));
    else
        %right branch of the solution
        Te(i)= TR + G2*erfc(xe(i)/2/sqrt(kR*time));
    end
end
hold on
plot(xe,Te,'r-') %plot exact solution
plot(s,Tc,'r+')  %plot the interface







