function [x_min, N] = lagrange(a0, b0, c0, epsilon, gamma, Nmax, fun)
%LAGRANGE  One-dimensional optimization using parabolic interpolation.
%   [x_min, N] = LAGRANGE(a0, b0, c0, epsilon, gamma, Nmax, fun)
%   finds the approximate minimum of a function fun(x) by fitting a
%   quadratic polynomial through three points and iteratively updating
%   the interval around the estimated minimum.
%
%   Inputs:
%     a0, b0, c0 - initial points (a0 < c0 < b0 recommended)
%     epsilon    - tolerance for the distance between a and b
%     gamma      - tolerance for change in d (approx. minimum)
%     Nmax       - maximum number of iterations
%     fun        - handle to the function to minimize
%
%   Outputs:
%     x_min      - estimated minimum of fun(x)
%     N          - number of iterations performed
%

    % Preallocate arrays for efficiency
    a = zeros(1, Nmax);
    b = zeros(1, Nmax);
    c = zeros(1, Nmax);
    d = zeros(1, Nmax);

    % Initialize points
    a(1) = a0;
    b(1) = b0;
    c(1) = c0;

    % Main iteration loop
    for i = 1:Nmax

        % Compute vertex (minimum) of the parabola through a, b, c
        d(i) = 0.5 * ( ...
            fun(a(i)) * (c(i)^2 - b(i)^2) + ...
            fun(c(i)) * (b(i)^2 - a(i)^2) + ...
            fun(b(i)) * (a(i)^2 - c(i)^2) ) / ...
            (fun(a(i)) * (c(i) - b(i)) + ...
             fun(c(i)) * (b(i) - a(i)) + ...
             fun(b(i)) * (a(i) - c(i)) );

        % Update points depending on where d lies
        if a(i) < d(i) && d(i) < c(i)
            if fun(d(i)) < fun(c(i))
                a(i+1) = a(i);
                b(i+1) = c(i);
                c(i+1) = d(i);
            else
                a(i+1) = d(i);
                b(i+1) = b(i);
                c(i+1) = c(i);
            end

        elseif c(i) < d(i) && d(i) < b(i)
            if fun(d(i)) < fun(c(i))
                a(i+1) = c(i);
                b(i+1) = b(i);
                c(i+1) = d(i);
            else
                a(i+1) = a(i);
                b(i+1) = d(i);
                c(i+1) = c(i);
            end

        else
            error('Algorithm is not convergent (points are inconsistent).');
        end

        % Stop if interval or function change is small enough
        if i > 1 && (abs(b(i+1) - a(i+1)) < epsilon || abs(d(i) - d(i-1)) < gamma)
            break
        end
    end

    % Final result
    x_min = d(i);
    N = i;
end