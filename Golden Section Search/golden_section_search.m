function [x_min, N] = golden_section_search(a0, b0, epsilon, N_max, fun)
%GOLDEN_SECTION_SEARCH  One-dimensional optimization using Golden Section.
%   [x_min, N] = GOLDEN_SECTION_SEARCH(a0, b0, epsilon, N_max, fun)
%   performs the Golden Section Search algorithm to locate the minimum
%   of a scalar function within the interval [a0, b0].
%
%   Inputs:
%     a0       - left endpoint of initial interval
%     b0       - right endpoint of initial interval
%     epsilon  - desired precision (stopping criterion)
%     N_max    - maximum number of iterations
%     fun      - handle to the function to minimize
%
%   Outputs:
%     x_min    - estimated position of the minimum
%     N        - number of iterations performed

    i = 0;
    alpha = (sqrt(5) - 1) / 2;  % golden ratio constant ≈ 0.618

    % Preallocate arrays for interval endpoints
    a = zeros(1, N_max);
    b = zeros(1, N_max);
    c = zeros(1, N_max);
    d = zeros(1, N_max);

    % Initialize first iteration
    a(1) = a0;
    b(1) = b0;
    c(1) = b(1) - alpha * (b(1) - a(1));
    d(1) = a(1) + alpha * (b(1) - a(1));

    % Main iteration loop
    while (b(i + 1) - a(i + 1)) >= epsilon
        % Compare function values at internal points
        if fun(c(i + 1)) < fun(d(i + 1))
            % Minimum lies in [a, d]
            a(i + 2) = a(i + 1);
            b(i + 2) = d(i + 1);
            d(i + 2) = c(i + 1);
            c(i + 2) = b(i + 2) - alpha * (b(i + 2) - a(i + 2));
        else
            % Minimum lies in [c, b]
            a(i + 2) = c(i + 1);
            b(i + 2) = b(i + 1);
            d(i + 2) = a(i + 2) + alpha * (b(i + 2) - a(i + 2));
            c(i + 2) = d(i + 1);
        end

        i = i + 1;

        % Safety check to avoid infinite loops
        if i > N_max
            error('Maximum number of iterations exceeded');
        end
    end

    % Approximate minimum is midpoint of final interval
    x_min = (a(i + 1) + b(i + 1)) / 2;
    N = i;
end
