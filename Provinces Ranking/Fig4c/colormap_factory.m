function C = colormap_factory(style, N)
% 生成 N 级离散同色系色带（接近你给的样例风格）。
% style: 'orange'|'purple'|'blue'|'green'|'red'|'yellow'|'teal'

if nargin<2, N=30; end
style = lower(string(style));
switch style
    case "orange", baseH = 28/360;  s0=0.25; s1=0.95; v0=1.00; v1=0.55;
    case "purple", baseH = 270/360; s0=0.20; s1=0.90; v0=0.98; v1=0.55;
    case "blue",   baseH = 210/360; s0=0.05; s1=0.85; v0=0.98; v1=0.55;
    case "green",  baseH = 130/360; s0=0.20; s1=0.90; v0=0.98; v1=0.55;
    case "red",    baseH =   5/360; s0=0.25; s1=0.95; v0=0.98; v1=0.55;
    case "yellow", baseH =  48/360; s0=0.05; s1=0.80; v0=1.00; v1=0.85;
    case "teal",   baseH = 170/360; s0=0.05; s1=0.85; v0=0.98; v1=0.60;
    otherwise, error('Unknown style: %s', style);
end
% 从浅到深（接近样例“上浅下深”）
s = linspace(s0, s1, N)'; v = linspace(v0, v1, N)';
h = baseH * ones(N,1);
C = hsv2rgb([h, s, v]);
end
