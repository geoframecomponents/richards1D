% BTCS scheme for the 2D Richards equation using the nested Newton method of
% Casulli & Zanolli.
clear all
close all
clc
global alpha thetas thetar n m Ks psic
global di dx dy dt K IMAX JMAX KL KR
%Phisical model parameters in SI units
day = 24*3600;
Ks = 0.062/day; %[meter/second]
thetas = 0.41;  %[-] saturated water content
thetar = 0.095; %[-] residuel water content
n = 1.31;       %[-] parameter n
m = 1-1/n;      %[-] parameter m
alpha = 1.9;    %[m^(-1)]
%critical value of psi where the maximum of the derivative is located:
psic = -1/alpha*((n-1)/n)^(1/n); 
%Domain
xL = 0;         %bottom
xR = 2;         %surface
yL = -1;        %left
yR = +1;        %right
IMAX = 50;      %number of control volumes
JMAX = 10;
dx = (xR-xL)/IMAX;      %mesh spacing
dy = (yR-yL)/JMAX;
x = linspace(xL+dx/2,xR-dx/2,IMAX);
y = linspace(yL+dy/2,yR-dy/2,JMAX);
tend = 3e5;             %set the final simulation time
time = 0;               %initial time
%set the initial condition (the right temperature-than it start freezing from the left)
for i=1:IMAX
    for j=1:JMAX
        psi(i,j) = -x(i);      % hydrostatic pressure  
    end
end
NMAX=100000;
for nit=1:NMAX
    dt = 1000;     %BTCS is unconditionalli stable (choose what you want)
    if(time+dt>tend)
        dt=tend-time;
    end
    %right (upper) boundary condition for pressure
    if(time<=1e5)
        psiR = -0.05 + 0.03*sin(2*pi*time/1e5);
    elseif(time>1e5 && time <= 1.8e5)
        psiR = +0.1;
    else
        psiR = -0.05 + 2952.45*exp(-time/18204.8);
    end
    %left (lower) boundary condition for pressure
    psiL = 0;
    KR = kappa(psiR);
    KL = kappa(psiL);
    %Compute lambda (heat conduction coefficient) and the internal energy at the old time level
    for i=1:IMAX
        for j=1:JMAX
            theta(i,j) = Thetaf(psi(i,j));
            K(i,j)     = kappa(psi(i,j));
        end
    end
    surf(x,y,psi')
    title(sprintf('Current time = %f',time))
    drawnow
    if(time>=tend)
        break
    end
    %compute right hand side and the linear part of the system
    for i=1:IMAX
        for j=1:JMAX
            % x (vertical) fluxes 
            if(i==1)
                Kp = 0.5*(K(i,j)+K(i+1,j));
                Km = 0.5*(K(i,j)+KL);
                %Dirichlet (pressure) boundary condition
                rhs(i,j) = theta(i,j) + dt/dx*(Kp-Km) + 2*Km*dt/dx^2*psiL;
            elseif(i==IMAX)
                Kp = 0.5*(K(i,j)+KR);
                Km = 0.5*(K(i,j)+K(i-1,j));
                %Dirichlet (pressure) boundary condition
                rhs(i,j)=theta(i,j) + dt/dx*(Kp-Km) + 2*Kp*dt/dx^2*psiR;
            else
                Kp = 0.5*(K(i,j)+K(i+1,j));
                Km = 0.5*(K(i,j)+K(i-1,j));
                rhs(i,j) = theta(i,j) + dt/dx*(Kp-Km);
            end
            % y fluxes with homogeneous Neumann BC
            % = do nothing!
        end
    end
    %now we start with the nested Newton method
    tol=1e-6;    %-12
    % initial guess
    psi=min(psi,psic);
    for iNewton=1:100 %outer Newton iterations
        %The task of the outer iterations is ONLY to linearize one of the
        %two nonlinear functions q1 or q2  (è il ciclo in k sul quaderno)
        di = zeros(IMAX,JMAX); %set the derivative of the nonlinear function to 0
        Mpsi = matop2D(psi);
        for i=1:IMAX
            for j=1:JMAX
                %compute the true non linear function f(psi)
                f(i,j)=Thetaf(psi(i,j))+ Mpsi(i,j)-rhs(i,j);
            end
        end
        outres = sqrt(sum(sum(f.*f))); %outer residual
        disp(sprintf('Outer iteration %d, outres=%e',iNewton, outres));
        if(outres<tol)
            break  %tolerance has been reached
        end
        psik = psi;          % save the value at the current outer iteration
        psi = max(psi,psic); % initial guess for the inner iteration
        for inner =1:100
            di = zeros(IMAX,JMAX);
            Mpsi = matop2D(psi);
            for i=1:IMAX
                for j=1:JMAX
                    fk(i,j)=Theta1(psi(i,j))-(Theta2(psik(i,j))+dTheta2(psik(i,j))*(psi(i,j)-psik(i,j)))+ Mpsi(i,j)-rhs(i,j);  %Tk is frozen
                    di(i,j)=dTheta1(psi(i,j))-dTheta2(psik(i,j));
                end
            end
            inres = sqrt(sum(sum(fk.*fk)));
            disp(sprintf(' -Inner iteration %d, inres= %e', inner,inres));
            if(inres<tol)
                break
            end
            dpsi = CGop2D(fk); % Can not use Thomas any more, because the matrix M is 5-diagonal!!
            psi = psi-dpsi;    %update temperature at the inner iteration
        end
        %in the outer iteration we freeze the temperature that we need for the
        %linearization of the second function
        
    end
    time = time+dt;     %advance time
end










