function [a, b] = contraction(xm1, xi, xp1, beta, N_max, fun)
%CONTRACTION  Adaptive bracket contraction around a local minimum.
%   [a, b] = CONTRACTION(xm1, xi, xp1, beta, N_max, fun)
%   contracts an existing bracket [xm1, xp1] around a presumed local minimum
%   near xi by shrinking the side where the function keeps decreasing.
%
%   Inputs:
%     xm1     - left endpoint of the current bracket
%     xi      - current middle point (near a local minimum)
%     xp1     - right endpoint of the current bracket
%     beta    - contraction factor in (0,1); steps scale like beta^j
%     N_max   - maximum number of contraction steps
%     fun     - handle to scalar objective function fun(x)
%
%   Outputs:
%     a, b    - updated bracket after contraction (a <= b)
%
%   Notes:
%     - Requires xm1 < xi < xp1.
%     - The algorithm decides contraction direction automatically.

    if beta <= 0 || beta >= 1
        error('beta must be in (0,1)');
    end
    if ~(xm1 < xi && xi < xp1)
        error('Require xm1 < xi < xp1');
    end
    if ~(isscalar(N_max) && N_max == floor(N_max) && N_max > 0)
        error('N_max must be a positive integer');
    end

    % Evaluate function values
    f_left  = fun(xm1);
    f_mid   = fun(xi);
    f_right = fun(xp1);

    if ~(f_mid <= f_left && f_mid <= f_right)
        warning('xi is not a local minimum relative to neighbors; results may be inaccurate.');
    end

    % Decide which side to contract
    if f_left < f_right
        % Contract the left side toward xi
        d = xi - xm1;
        j = 1;
        x_new = xi - beta^j * d;
        f_new = fun(x_new);
        while f_new < f_mid
            xm1 = x_new;
            j = j + 1;
            if j > N_max, break; end
            x_new = xi - beta^j * d;
            f_new = fun(x_new);
        end
    else
        % Contract the right side toward xi
        d = xp1 - xi;
        j = 1;
        x_new = xi + beta^j * d;
        f_new = fun(x_new);
        while f_new < f_mid
            xp1 = x_new;
            j = j + 1;
            if j > N_max, break; end
            x_new = xi + beta^j * d;
            f_new = fun(x_new);
        end
    end

    % Return ordered interval
    a = min(xm1, xp1);
    b = max(xm1, xp1);
end
