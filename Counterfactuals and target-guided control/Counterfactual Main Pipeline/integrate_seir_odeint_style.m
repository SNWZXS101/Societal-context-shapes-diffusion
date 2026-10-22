function out = integrate_seir_odeint_style(f, init, beta, gamma)
% t = 0:1:T-1；f 采用“最近邻（ties-to-even）”取样
% 稳定性策略：rhs 里以非负裁剪参与 inc=βSI/N 的计算；可选 NonNegative

T = numel(f);
tspan = 0:(T-1);
y0 = [max(0,init.Se0); max(0,init.S0); max(0,init.I0); max(0,init.R0)];

useNonneg = true;                       % 如需完全复刻最初口径，可设为 false

    function idx = idx_from_t(tt)
        rr  = round_half_to_even(tt);
        rr  = max(0, min(T-1, rr));
        idx = rr + 1;
    end
    function val = NN_f(tt), val = f(idx_from_t(tt)); end

    function dydt = rhs(tt, y)
        % —— 用“非负裁剪”的状态参与传播项，避免负值污染右端 ——
        Se = max(0.0, y(1)); S = max(0.0, y(2));
        I  = max(0.0, y(3)); R = max(0.0, y(4));
        den = S + I + R; if den <= 0, den = eps; end
        inc = beta * S * I / den;
        ft  = NN_f(tt);                  % 不再对 ft 做额外限幅
        dydt = [-ft; ft - inc; inc - gamma*I; gamma*I];
    end

% —— 关键变化：去掉 MaxStep/InitialStep，放宽容差，让自适应更自由 ——
opts = odeset('RelTol',1e-7,'AbsTol',1e-10);
if useNonneg, opts = odeset(opts,'NonNegative',1:4); end

[tt, YY] = ode45(@rhs, tspan, y0, opts);   % 仍输出整日网格

Se = max(0, YY(:,1)); S = max(0, YY(:,2));
I  = max(0, YY(:,3)); R = max(0, YY(:,4));
den = S + I + R; den(den<=0) = eps;
I_ode = beta .* S .* I ./ den;

out = struct('t',tt, 'Se',Se, 'S',S, 'I',I, 'R',R, 'I_ode',I_ode);
end

function r = round_half_to_even(x)
fx=floor(x); frac=x-fx; r=fx;
gt = frac>0.5; eq = abs(frac-0.5)<=1e-12;
r(gt)=fx(gt)+1; odd = mod(fx,2)~=0; bump = eq & odd; r(bump)=fx(bump)+1;
end
