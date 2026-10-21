function C = nice_colors(K)
base = [ ...
    0.1216 0.4667 0.7059;  % blue
    1.0000 0.4980 0.0549;  % orange
    0.1725 0.6275 0.1725;  % green
    0.8392 0.1529 0.1569;  % red
    0.5804 0.4039 0.7412;  % purple
    0.5490 0.3373 0.2941;  % brown
    0.8902 0.4667 0.7608;  % pink
    0.4980 0.4980 0.4980]; % gray
C = repmat(base, ceil(K/size(base,1)), 1);
C = C(1:K,:);
end
