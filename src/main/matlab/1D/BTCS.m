% BTCS scheme for the 1D Richards equation using the nested Newton method of
% Casulli & Zanolli.

clear all 
close all 
clc
global alpha thetas thetar n m Ks psic
%epsilon= 0.1; % regularization parameter: size of the regularization interval
%Phisical model parameters in SI units
day     = 24*3600;
Ks      = 0.062/day;    %[meter/second]
thetas  = 0.41;         %[-] saturated water content
thetar  = 0.095;        %[-] residuel water content
n       = 1.31;         %[-] parameter n
m       = 1-1/n;        %[-] parameter m
alpha   = 1.9;          %[m^(-1)]
psic    = -1/alpha*((n-1)/n)^(1/n);  %critical value of psi where the maximum of the derivative is located
%Domain
xL = 0;                 %bottom
xR = 2;                 %surface
IMAX = 100;             %number of control volumes
dx = (xR-xL)/IMAX;      %mesh spacing
x = linspace(xL+dx/2,xR-dx/2,IMAX);
tend = 3e8;             %set the final simulation time
time = 0;               %initial time
%set the initial condition (the right temperature-than it start freezing from the left)
for i=1:IMAX
    psi(i) = -x(i);      % hydrostatic pressure   
end
NMAX=100000;
for nit=1:NMAX
    dt = 100000;   %BTCS is unconditionalli stable (choose what you want)
    if(time+dt>tend)
        dt=tend-time;
    end
    %right boundary condition for pressure
    if(time<=1e5)
        psiR = -0.05+0.03*sin(2*pi*time/1e5);
    elseif(time>1e5 && time<=1.8e5)
        psiR = +0.1;
    else
        psiR = -0.05+2952.45*exp(-time/18204.8);
    end
    %left
    psiL = 0;
    KR = kappa(psiR); 
    KL = kappa(psiL);
    %Compute lambda (heat conduction coefficient) and the internal energy at the old time level
    for i=1:IMAX
        theta(i)=Thetaf(psi(i));
        K(i) = kappa(psi(i));           
    end
    subplot(2,1,1)
    plot(x,psi,'o') 
    subplot(2,1,2)
    plot(x,theta,'o')      
    title(sprintf('Current time = %f',time))
    drawnow
    if(time>=tend)
        break
    end    
    %compute right hand side and the linear part of the system
    
    for i=1:IMAX
        if(i==1)
            Kp = 0.5*(K(i)+K(i+1));
            Km = 0.5*(K(i)+KL);
            %dirichlet (pressure) boundary condition
            a(i) = 0;
            b(i) = dt/dx^2*(2*Km+Kp);
            c(i) = -Kp*dt/dx^2;
            rhs(i) = theta(i) + dt/dx*(Kp-Km) + 2*Km*dt/dx^2*psiL; 
        elseif(i==IMAX)
            Kp = 0.5*(K(i)+KR);
            Km = 0.5*(K(i)+K(i-1));
            %dirichlet (pressure boundary condition
            a(i) = -Km*dt/dx^2;
            b(i) = dt/dx^2*(Km+2*Kp);
            c(i) = 0;
            rhs(i) = theta(i) + dt/dx*(Kp-Km) + 2*Kp*dt/dx^2*psiR;              
        else
            Kp = 0.5*(K(i)+K(i+1));
            Km = 0.5*(K(i)+K(i-1));
            a(i) = -Km*dt/dx^2;
            b(i) = dt/dx^2*(Km+Kp);
            c(i) = -Kp*dt/dx^2;
            rhs(i) = theta(i) + dt/dx*(Kp-Km);           
        end
    end
    %now we start with the nested Newton method
    tol = 1e-12;
    % initial guess
    psi = min(psi,psic);
    for iNewton=1:100 %outer Newton iterations
        %The task of the outer iterations is ONLY to linearize one of the
        %two nonlinear functions q1 or q2  (� il ciclo in k sul quaderno)
        for i=1:IMAX
            %compute the true non linear function f(psi)
            f(i) = Thetaf(psi(i))-rhs(i);        
            if(i==1)
                f(i) = f(i) + b(i)*psi(i) + c(i)*psi(i+1);
            elseif(i==IMAX)
                f(i) = f(i)+a(i)*psi(i-1)+b(i)*psi(i);
            else
                f(i) = f(i)+a(i)*psi(i-1)+b(i)*psi(i)+c(i)*psi(i+1);
            end
        end
        outres = sqrt(sum(f.*f)); %outer residual
        disp(sprintf('Outer iteration %d, outres=%e',iNewton, outres));    
        if(outres<tol)
            break %tolerance has been reached
    end

    psik = psi;             % save the value at the current outer iteration    
    psi = max(psi,psic); % initial guess for the inner iteration

    for inner=1:100
        for i=1:IMAX
            fk(i)=Theta1(psi(i))-(Theta2(psik(i))+dTheta2(psik(i))*(psi(i)-psik(i)))-rhs(i);  %Tk is frozen
            if(i==1)
                fk(i)=fk(i)+b(i)*psi(i)+c(i)*psi(i+1);
            elseif(i==IMAX)
                fk(i)=fk(i)+a(i)*psi(i-1)+b(i)*psi(i);
            else
                fk(i)=fk(i)+a(i)*psi(i-1)+b(i)*psi(i)+c(i)*psi(i+1);
            end
            di(i)=dTheta1(psi(i))-dTheta2(psik(i));
        end
        inres=sqrt(sum(fk.*fk));
        disp(sprintf(' -Inner iteration %d, inres= %e', inner,inres));
        if(inres<tol)
            break
        end
        dpsi = Thomas(a,b+di,c,fk); %inner Newton step
        psi = psi(:)-dpsi(:);
    end    
end
    time = time+dt;     %advance time
end










