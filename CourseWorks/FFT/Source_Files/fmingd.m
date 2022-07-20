function sol = fmingd(f,x0,grad)
%FMINGD finds local minimum of unconstrained function f via gradient decent
%method.
%   Detailed explanation goes here
%f - function handle of optimized function, 
%grad - function handle of gradient of f
%x0 - initial point
%sol - output structure: sol.xmin - minimized x, sol.fmin - local min of f
%                        sol.stps - vector of steps coordinates 
%                        sol.stps_val - f values on each step

    alpha = 0.1; rel_tol = 0.001;
    xk1 = x0; fk1 = f(x0);
    sol.xmin = xk1; sol.fmin = fk1;
    sol.stps = xk1; sol.stps_val = fk1;
    
    while true      %do .. while imitation
        xk = xk1; fk = fk1;
        xk1 = rec_f(xk); fk1 = f(xk1); %x_k+1 = rec(x_k) ...
        sol.stps = cat(2,sol.stps,xk1);
        sol.stps_val = cat(2,sol.stps_val,fk1);
        if (abs(f(xk1)-f(xk)) < rel_tol) 
            sol.xmin = xk1; sol.fmin = fk1;
            break
        end
    end

function xk1 = rec_f(xk)
    xk1 = xk - alpha*(grad(xk))';
end
end

