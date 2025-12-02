function [x_min, N] = powell(x0, d, eps, Nmax, f)
% Powell's method
% using fminsearch
%
% f    - objective function handle
% x0   - initial point (column vector)
% d    - initial directions matrix (each column is a direction d_j)
% eps  - tolerance
% Nmax - maximum number of iterations
%
% Returns:
%   x_min - approximate minimizer
%   N     - number of outer iterations performed

    n = size(d, 2);
    x = x0;

    for i = 1:Nmax
        P = zeros(length(x), n+1);
        P(:,1) = x;
        
        for j = 1:n
            % Line search
            phi = @(a) f(P(:,j) + a*d(:,j));
            
            alpha = fminsearch(phi, 0);
            
            P(:,j+1) = P(:,j) + alpha * d(:,j);
        end
        
        p_n = P(:,n+1);
        
        % Stopping condition
        if abs(f(p_n) - f(x)) < eps
            x_min = p_n;
            N = i;
            return
        end
        
        % Update directions
        for j = 1:n-1
            d(:,j) = d(:,j+1);
        end
        
        d(:,n) = p_n - P(:,1);
        
        phi_new = @(a) f(p_n + a*d(:,n));
        alpha_new = fminsearch(phi_new, 0);
        
        p_np1 = p_n + alpha_new * d(:,n);
        
        x = p_np1;
    end

    % Return error
    error("Maximum number of iterations exceeded");
end
