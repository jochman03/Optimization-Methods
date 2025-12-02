function [x_min, N] = penalty_method(x0, c1, c, eps, Nmax, f, s)
%PENALTY_METHOD.
%   [x_min, N] = PENALTY_METHOD(x0, c1, c, eps, Nmax, f, s)
%   minimizes an objective function f(x) subject to constraints encoded
%   by the penalty function s(x).
%
%   Inputs:
%     x0   - initial point from which the algorithm starts
%     c1   - initial value of the penalty coefficient c_i
%     c    - update function for the penalty coefficient
%     eps  - stopping tolerance for changes in successive iterates
%     Nmax - maximum allowed number of iterations
%     f    - handle to the objective function f(x)
%     s    - handle to the penalty function s(x),
%
%   Outputs:
%     x_min - estimated minimum of the constrained optimization problem
%     N     - number of outer iterations performed
%
%   Notes:
%     - The algorithm solves a sequence of unconstrained problems:
%           minimize  F_i(x) = f(x) + c_i * s(x)
%       where c_i changes at each iteration.
%
%     - Ideally, s(x) penalizes constraint violation, e.g.:
%           s(x) = sum( max(0, g_j(x))^2 )   for inequality constraints
%
%     - Each unconstrained subproblem is solved using fminunc().
%
%     - The method stops when:
%           || x_{i} - x_{i-1} || < eps
%       or when the number of iterations reaches Nmax.
%

    ci = c1;
    x_last = x0;

    for i = 1:Nmax
        % penalized objective
        F = @(x) f(x) + ci * s(x);
        x_star = fminunc(F, x_last, optimoptions('fminunc','Display','off'));
        
        % stopping condition
        if norm(x_star - x_last, 2) < eps
            x_min = x_star;
            N = i;
            return;
        end

        x_last = x_star;
        % update penalty coefficient
        ci = c(ci);  
    end

    error("Maximum number of iterations exceeded");
end