function [a, b] = expansion(x0, x1, alpha, N_max, fun)
%   Bracketing expansion starting from an arbitrary x0.
%   [a, b] = expansion(x0, x1, alpha, N_max, fun)
%   reduces the problem to y = x - x0, performs geometric expansion
%   around 0 to find an interval [a,b] such that fun has a local
%   minimum inside (standard bracketing step), then translates back.
%
%   Inputs:
%     x0       - reference point (start)
%     x1       - initial trial point (x1 ~= x0 recommended)
%     alpha    - expansion factor (> 1)
%     N_max    - maximum number of expansion iterations
%     fun      - handle to scalar objective function fun(x)
%
%   Outputs:
%     a, b     - endpoints of the bracket in the original x-domain

    % Sanity check
    if alpha <= 1
        error('alpha must be > 1');
    end

    % Reduce to y-domain: y = x - x0, so we expand around 0
    g  = @(y) fun(y + x0);
    y0 = 0;
    y1 = x1 - x0;

    % Evaluate objective at start and first trial
    f0 = g(y0);
    f1 = g(y1);

    % If values are equal at the start, the trivial bracket is [min, max]
    if f0 == f1
        a = min(x0, x1);
        b = max(x0, x1);
        return
    end

    % If we moved uphill, flip direction (use symmetry around 0)
    if f1 > f0
        y1 = -y1;
        f1 = g(y1);
        % If still not downhill, the symmetric points bracket the minimum
        if f1 >= f0
            ay = min(y1, -y1);
            by = max(y1, -y1);
            a = ay + x0;
            b = by + x0;
            return
        end
    end

    % Geometric expansion from 0 in the current downhill direction
    i     = 1;
    yi_1  = y0;  fi_1  = f0;         % x^{(i-1)}
    yi    = y1;  fi    = f1;         % x^{(i)}
    yip1  = (alpha^i) * y1;          % x^{(i+1)} = alpha^i * y1
    fip1  = g(yip1);

    % Keep expanding while the function keeps decreasing
    while fi > fip1
        if i >= N_max
            error('Maximum number of iterations exceeded (N_max).');
        end
        yi_1 = yi;   fi_1 = fi;
        yi   = yip1; fi   = fip1;

        i    = i + 1;
        yip1 = (alpha^i) * y1;
        fip1 = g(yip1);
    end

    % Translate the final y-bracket back to x-domain
    ay = min(yi_1, yip1);
    by = max(yi_1, yip1);
    a  = ay + x0;
    b  = by + x0;
end
