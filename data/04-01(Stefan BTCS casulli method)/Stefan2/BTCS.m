% BTCS scheme for the Stefan problem using the nested Newton method of
% Casulli & Zanolli.

clear all 
close all 
clc
global kL kR hc cL cR Tc TR TL rhoL rhoR lambdaL lambdaR epsilon
epsilon= 0.1; % regularization parameter: size of the regularization interval
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
   
end
NMAX=100000;

for n=1:NMAX
    dt = 3600;
    if(time+dt>tend)
        dt=tend-time;
    end
    
    %compute lambda (heat conduction coefficient) and the internal energy
    % at the old time level
    for i=1:IMAX
        qn(i)=q(T(i));
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
    
    %compute right hand side and the linear part of the system
    for i=1:IMAX
        if(i==1)
            lambdap = 0.5*(lambda(i)+lambda(i+1));
            lambdam = lambda(i);
            a(i)=0; %lower diagonal of matrix M
            b(i)=dt/dx^2*(2*lambdam+lambdap);%diagonal
            c(i)= -lambdap*dt/dx^2;%upper diagonal
            rhs(i)= qn(i)+2*lambdam*dt/dx^2*TL;
        elseif(i==IMAX)
            lambdap = lambda(i);
            lambdam = 0.5*(lambda(i)+lambda(i-1));
            a(i)=-lambdam*dt/dx^2;
            b(i)=dt/dx^2*(lambdam+2*lambdap);
            c(i)=0;
            rhs(i)=qn(i)+2*lambdap*dt/dx^2*TR;
            
           
        else
            lambdap = 0.5*(lambda(i)+lambda(i+1));
            lambdam = 0.5*(lambda(i)+lambda(i-1));
            a(i)=-lambdam*dt/dx^2;
            b(i)=dt/dx^2*(lambdam+lambdap);
            c(i)=-lambdap*dt/dx^2;
            rhs(i)=qn(i);
           
        end
    end
    %now we start with the nested Newton method
    tol=1e-12*rhoR*hc; %rescaled tolerance(round -off errors!!)
    % initial guess
    T=min(T,Tc-epsilon);
    for iNewton=1:100 %outer Newton iterations
        %The task of the outer iterations is ONLY to linearize one of the
        %two nonlinear functions q1 or q2  (è il ciclo in k sul quaderno)
    for i=1:IMAX
        %compute the true non linear function f(T)
        f(i)=q(T(i))-rhs(i);
        
        if(i==1)
            f(i)=f(i)+b(i)*T(i)+c(i)*T(i+1);
        elseif(i==IMAX)
            f(i)=f(i)+a(i)*T(i-1)+b(i)*T(i);
        else
            f(i)=f(i)+a(i)*T(i-1)+b(i)*T(i)+c(i)*T(i+1);
        end
    end
    outres=sqrt(sum(f.*f)); %outer residual
    disp(sprintf('Outer iteration %d, outres=%e',iNewton, outres));
    
    if(outres<tol)
        break %tolerance has been reached
    end
    Tk=T; % save the value at the current outer iteration
    
    T=max(Tk,Tc-epsilon); % initial guess for the inner iteration
    for inner =1:100
       for i=1:IMAX
          fk(i)=q1(T(i))-(q2(Tk(i))+dq2(Tk(i))*(T(i)-Tk(i)))-rhs(i);  %Tk is frozen 
           if(i==1)
            fk(i)=fk(i)+b(i)*T(i)+c(i)*T(i+1);
        elseif(i==IMAX)
            fk(i)=fk(i)+a(i)*T(i-1)+b(i)*T(i);
        else
            fk(i)=fk(i)+a(i)*T(i-1)+b(i)*T(i)+c(i)*T(i+1);
        end
        di(i)=dq1(T(i))-dq2(Tk(i)); %dq1 non decreasing and dq2(Tk(i)) is constant
       end
        inres=sqrt(sum(fk.*fk));
        disp(sprintf(' -Inner iteration %d, inres= %e', inner,inres));
        if(inres<tol)
            break
        end
        dT = Thomas(a,b+di,c,fk); %inner Newton step
        T=T(:)-dT(:); %update temperature at the inner iteration
      end
   %in the outer iteration we freeze the temperature that we need for the
   %linearization of the second function 
    
end
    time = time+dt;     %advance time
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







