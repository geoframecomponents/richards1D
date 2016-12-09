function Apsi = matop2D(psi)
global di dx dy dt K IMAX JMAX KL KR
% diagonal part including the derivatives of the nonlinear functions
Apsi = di.*psi;
% linear part
for i=1:IMAX
    for j=1:JMAX
        % x fluxes
        if(i==1)
            Kp = 0.5*(K(i,j)+K(i+1,j));
            Km = 0.5*(K(i,j)+KL);
            Apsi(i,j) = Apsi(i,j)-dt/dx^2*( Kp*(psi(i+1,j)-psi(i,j))...
                                           - 2*Km*(psi(i,j)-0) );
        elseif(i==IMAX)
            Kp = 0.5*(K(i,j) + KR);
            Km = 0.5*(K(i,j) + K(i-1,j));
            Apsi(i,j) = Apsi(i,j)-dt/dx^2*( 2*Kp*(0-psi(i,j))...
                                           -Km*(psi(i,j)-psi(i-1,j)) );
        else
            Kp = 0.5*(K(i,j) + K(i+1,j));
            Km = 0.5*(K(i,j) + K(i-1,j));
            Apsi(i,j) = Apsi(i,j)-dt/dx^2*( Kp*(psi(i+1,j)-psi(i,j))...
                                           -Km*(psi(i,j)-psi(i-1,j)) );
        end
        % y fluxes
        if (j==1)
           Kp = 0.5*(K(i,j) + K(i,j+1));
           Apsi(i,j) = Apsi(i,j) - dt/dy^2*( Kp*(psi(i,j+1)-psi(i,j)) - 0 );
        elseif (j==JMAX)
           Km = 0.5*(K(i,j)+K(i,j-1));
           Apsi(i,j) = Apsi(i,j) - dt/dy^2*( 0 - Km*(psi(i,j)-psi(i,j-1)) );
        else
           Kp = 0.5*(K(i,j)+K(i,j+1));
           Km = 0.5*(K(i,j)+K(i,j-1));
           Apsi(i,j) = Apsi(i,j) - dt/dy^2*( Kp*(psi(i,j+1)-psi(i,j))...
                                            -Km*(psi(i,j)-psi(i,j-1)) );
        end
    end
end
