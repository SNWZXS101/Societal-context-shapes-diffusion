
function Ntot = try_get_N(W, i)
Ntot = [];
if ~isempty(W.init)
    ini = W.init{i};
    if isstruct(ini)
        have = isfield(ini,'S') && isfield(ini,'I') && isfield(ini,'R');
        if have
            S0 = double(ini.S(1)); I0 = double(ini.I(1)); R0 = double(ini.R(1));
            Ntot = S0 + I0 + R0;
        end
    elseif isnumeric(ini) && numel(ini)>=3
        Ntot = sum(ini(1:3));
    end
end
end