
function s = nansum(x), x=x(:); x(isnan(x))=[]; s=sum(x); end