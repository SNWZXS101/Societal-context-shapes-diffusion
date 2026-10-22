function out = integrate_seir_odeint_style(f, init, beta, gamma)
% INTEGRATE_SEIR_ODEINT_STYLE
% Replicates the Python odeint setup in your code:
% - Time grid: t = 0,1,2,...,T-1 (days)
% - External input: NN_f(t) = f( round_half_to_even(t), clamped to [0,T-1] )
% Returns a struct with fields: t, Se, S, I, R, I_ode

T = numel(f);
tspan = 0:(T-1);
y0 = [init.Se0; init.S0; init.I0; init.R0];

    function idx = idx_from_t(tt)
        % Python's round (banker's rounding) with clamping, then 1-based
        rr  = round_half_to_even(tt);
        rr  = max(0, min(T-1, rr));
        idx = rr + 1;
    end

    function val = NN_f(tt)
        val = f(idx_from_t(tt));
    end

    function dydt = rhs(tt, y)
        Se = y(1); S = y(2); I = y(3); R = y(4);
        denom = S + I + R;
        if denom <= 0, denom = eps; end
        inc = beta * S * I / denom;
        ft  = NN_f(tt);
        dydt = [ -ft;
                  ft - inc;
                  inc - gamma * I;
                  gamma * I ];
    end

% 精细容差 + 小步长限制，尽量贴近 LSODA 的行为，同时保证整数日输出
opts = odeset('RelTol',1e-9,'AbsTol',1e-12,'MaxStep',0.1,'InitialStep',0.1);
[tt, YY] = ode45(@rhs, tspan, y0, opts);

Se = YY(:,1); S = YY(:,2); I = YY(:,3); R = YY(:,4);
den = S + I + R; den(den<=0) = eps;
I_ode = beta .* S .* I ./ den;

out = struct('t',tt, 'Se',Se, 'S',S, 'I',I, 'R',R, 'I_ode',I_ode);
end

% -------- helpers --------
function r = round_half_to_even(x)
% Banker's rounding for real vectors (ties-to-even).
fx   = floor(x);
frac = x - fx;
r    = fx;
gt   = frac > 0.5;
eq   = abs(frac - 0.5) <= 1e-12;    % tolerance on halves
r(gt) = fx(gt) + 1;
odd   = mod(fx,2) ~= 0;
bump  = eq & odd;                    % 0.5 above odd -> ceil
r(bump) = fx(bump) + 1;
% eq & even -> stay at fx
end
