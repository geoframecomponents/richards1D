%Newton method for the solution of the nonlinear algebric equation g(gamma)=0

function gamma = Newton(gamma0)
gamma = gamma0; %initial guess
tol = 1e-12;    %tollerance (if res<tol can stop
for iNewton = 1:100
    gk = g(gamma);             %function of the current N iteration
    res = abs(gk);             %residual
    disp(sprintf('- Iteration %d, res= %e',iNewton,res))
    if(res<tol)         %if res is below tol
        break           %than stop the iteration
    end
    dgamma= -gk/dg(gamma);    %compute the Newton step
    gamma= gamma + dgamma;    %solution at the next iteration
end