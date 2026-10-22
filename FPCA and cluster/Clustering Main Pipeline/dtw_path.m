function [p,q] = dtw_path(x, y)
% return warping path indices (1-based)
x = x(:); y = y(:);
nx=numel(x); ny=numel(y);
D = (x - y').^2;
C = inf(nx+1,ny+1); C(1,1)=0;
phi = zeros(nx,ny); % 1:diag, 2:up, 3:left
for i=1:nx
    for j=1:ny
        [m,ix] = min([C(i,j), C(i, j+1), C(i+1, j)]);
        C(i+1,j+1) = D(i,j) + m;
        phi(i,j) = ix;
    end
end
% backtrack
i=nx; j=ny; P=[]; Q=[];
while i>0 && j>0
    P=[i;P]; Q=[j;Q]; %#ok<AGROW>
    if phi(i,j)==1 
        i=i-1; j=j-1;
    elseif phi(i,j)==2
        i=i-1;
    else
        j=j-1; 
    end
end
p=P; q=Q;
end
