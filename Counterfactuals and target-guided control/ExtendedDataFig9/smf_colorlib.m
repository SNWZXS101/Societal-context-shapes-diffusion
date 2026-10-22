function C = smf_colorlib(which, varargin)
% SMF_COLORLIB  Nature级统一配色库（分类/顺序/发散）+ 可选变暗/增亮/去饱和
% 用法：
%   C = smf_colorlib('groups8');                 % 返回预定义调色板（Nx3, 0-1）
%   C = smf_colorlib('groups8','darken',0.15);   % 整体加深
%   C = smf_colorlib('viridis', 'K', 11);        % 从连续色带抽取K=11个离散色
%   C = smf_colorlib('categories30', 'K', 30);   % 生成30个互异分类色（算法生成）
%
% 也提供简便函数：
%   Ck = smf_color_pick('okabeito8', 12);        % 超出内置长度时自动扩展生成

% ------------ 解析参数 ------------
p = inputParser;
p.addParameter('darken', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<1);
p.addParameter('lighten', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<1);
p.addParameter('desaturate', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('K', [], @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.parse(varargin{:});
dval = p.Results.darken;
lval = p.Results.lighten;
sval = p.Results.desaturate;
K    = p.Results.K;

if nargin==0 || isempty(which), which='groups8'; end
key = lower(string(which));

% 别名处理（_deep = 默认加深）
DEFAULT_DEEP = 0.15;
if endsWith(key,'_deep')
    key = extractBefore(key, strlength(key)-4);
    dval = max(dval, DEFAULT_DEEP);
end

% ------------ 主开关：返回颜色矩阵或连续色带采样 ------------
switch key


    case 'RZJS7'
        C = [ ...
            0.961 0.773 0.094;   % 金黄  (#F5C518)
            0.961 0.624 0.000;   % 番红  (#F59F00)
            0.906 0.435 0.318;   % 珊瑚  (#E76F51)
            0.784 0.333 0.541;   % 玫红  (#C8558A)
            0.424 0.361 0.906;   % 紫罗  (#6C5CE7)
            0.247 0.471 0.710;   % 群青  (#3F78B5)
            0.184 0.435 0.494];  % 灰蓝  (#2F6F7E)


        % ===== 主打配色（来自两张图片）=====
    case 'maina7'
        C = [ ...
            79 177  178;
            51 123  172;
            144 180 207;
            169 217 187;
            202 235 231;
            171 215 223;
            213 232 241] / 255;

    case 'maina8'   % 图1：152,216,59 → 70,46,122（青绿到靛紫，8色）
        C = [ ...
            152 216  59;
            96 201  96;
            48 180 123;
            28 153 135;
            37 129 142;
            47 105 142;
            61  77 137;
            70  46 122] / 255;

    case 'mainb8'   % 图2：取前8色（绿→黄→橙→红，8色）
        C = [ ...
            0  99  53;
            0 157  83;
            119 181 106;
            179 206 126;
            241 241 183;
            242 203 123;
            232 136  84;
            227  58  40] / 255;

    case 'maina9'   % 图2完整9色版本（含最深红 167,25,39）
        C = [ ...
            0  99  53;
            0 157  83;
            119 181 106;
            179 206 126;
            241 241 183;
            242 203 123;
            232 136  84;
            227  58  40;
            167  25  39] / 255;
        % ===== 扩展版：mainA30 / mainB30（亮度交替更明显）=====
    case 'maina30'
        base = smf_colorlib('mainA8');     % 8 基色（已在库中）
        C    = expand_palette_strong(base, 30);
    case 'mainb30'
        base = smf_colorlib('mainB9');     % 9 基色（已在库中）
        C    = expand_palette_strong(base, 30);
    

        % ===== 你原来的几套（保留）=====
    case 'clusters4'
        C = [0.12 0.47 0.71; 0.17 0.63 0.17; 0.58 0.40 0.74; 0.84 0.15 0.16];
    case 'clusters4_alt1'
        C = [0.11 0.62 0.52; 0.85 0.37 0.01; 0.46 0.44 0.70; 0.91 0.16 0.54];
    case 'clusters4_alt2'
        C = [0.10 0.28 0.46; 0.78 0.64 0.00; 0.13 0.55 0.28; 0.75 0.11 0.20];
    case 'clusters4_alt3'
        C = [0.00 0.59 0.53; 0.80 0.47 0.13; 0.48 0.31 0.64; 0.00 0.45 0.45];
    case 'variants3'
        C = [0.12 0.47 0.71; 1.00 0.50 0.05; 0.17 0.63 0.17];
    case 'waves'
        C = [0.23 0.48 0.60; 0.44 0.69 0.81; 0.64 0.86 0.95; 0.96 0.97 0.98];
    case 'marsh'
        C = [0.35 0.68 0.87; 0.56 0.81 0.94; 0.00 0.50 0.94; 0.63 0.85 0.84];
    case 'moss'
        C = [0.28 0.28 0.47; 0.24 0.54 0.65; 0.27 0.75 0.70; 0.67 0.82 0.35];
    case 'icecream'
        C = [0.67 0.39 0.47; 0.87 0.87 0.86; 0.99 0.86 0.84; 0.93 0.66 0.65];
    case 'sunset'
        C = [0.80 0.29 0.54; 0.96 0.38 0.65; 0.99 0.59 0.22; 0.93 0.67 0.47];
    case 'velvet'
        C = [0.07 0.20 0.40; 0.35 0.00 0.60; 0.59 0.00 0.59; 0.60 0.78 1.00];
    case 'soviet'
        C = [0.74 0.82 0.86; 0.98 0.98 0.98; 0.87 0.94 0.98; 0.27 0.32 0.45];
    case 'coconut'
        C = [0.04 0.55 0.44; 0.04 0.73 0.53; 0.58 0.85 0.66; 0.96 0.85 0.67];
    case 'root'
        C = [0.33 0.55 0.39; 0.99 0.79 0.52; 0.76 0.52 0.22; 0.92 0.70 0.89];

        % ===== 经典/常用（Nature友好）=====
    case 'classic10' % Tableau-like
        C = [ ...
            0.1216 0.4667 0.7059; ...
            1.0000 0.4980 0.0549; ...
            0.1725 0.6275 0.1725; ...
            0.8392 0.1529 0.1569; ...
            0.5804 0.4039 0.7412; ...
            0.5490 0.3373 0.2941; ...
            0.8902 0.4667 0.7608; ...
            0.4980 0.4980 0.4980; ...
            0.7373 0.7412 0.1333; ...
            0.0902 0.7451 0.8118];
    case 'okabeito8'  % 色盲友好（黑+7色）
        C = [ ...
            0.000 0.000 0.000;  % black
            0.902 0.624 0.000;  % orange (230,159,0)
            0.337 0.706 0.914;  % sky blue (86,180,233)
            0.000 0.620 0.451;  % bluish green (0,158,115)
            0.941 0.894 0.259;  % yellow (240,228,66)
            0.000 0.447 0.698;  % blue (0,114,178)
            0.835 0.369 0.000;  % vermillion (213,94,0)
            0.800 0.475 0.655]; % reddish purple (204,121,167)
    case 'groups8' % 我们项目的8组：Population/Healthcare/Economy/Urbanicity/Mobility/Tourism/Region/Climate
        C = [ ...
            0.1216 0.4667 0.7059; ... % 蓝 (Population)
            0.1725 0.6275 0.1725; ... % 绿 (Healthcare)
            0.8392 0.1529 0.1569; ... % 红 (Economy)
            0.5804 0.4039 0.7412; ... % 紫 (Urbanicity)
            1.0000 0.4980 0.0549; ... % 橙 (Mobility)
            0.0902 0.7451 0.8118; ... % 青 (Tourism)
            0.7373 0.7412 0.1333; ... % 橄榄 (Region)
            0.5490 0.3373 0.2941];    % 棕 (Climate)
    case 'gray'
        C = [0.5 0.5 0.5];

        % ===== ColorBrewer（质性）=====
    case 'cb_paired12'
        C = [ ...
            0.65098 0.80784 0.89020;
            0.12157 0.47059 0.70588;
            0.69804 0.87451 0.54118;
            0.20000 0.62745 0.17255;
            0.98431 0.60392 0.60000;
            0.89020 0.10196 0.10980;
            0.99216 0.74902 0.43529;
            1.00000 0.49804 0.00000;
            0.79216 0.69804 0.83922;
            0.41569 0.23922 0.60392;
            1.00000 1.00000 0.60000;
            0.69412 0.34902 0.15686];
    case 'cb_set3_12'
        C = [ ...
            0.55294 0.82745 0.78039;
            1.00000 1.00000 0.70196;
            0.74510 0.72941 0.85490;
            0.98431 0.50196 0.44706;
            0.50196 0.69412 0.82745;
            0.99216 0.70588 0.38431;
            0.70196 0.87059 0.41176;
            0.98824 0.80392 0.89804;
            0.85098 0.85098 0.85098;
            0.73725 0.50196 0.74118;
            0.80000 0.92157 0.77255;
            1.00000 0.92941 0.43529];
    case 'cb_dark2_8'
        C = [ ...
            0.10588 0.61961 0.46667;
            0.85098 0.37255 0.00784;
            0.45882 0.43922 0.70196;
            0.90588 0.16078 0.54118;
            0.40000 0.65098 0.11765;
            0.90196 0.67059 0.00784;
            0.65098 0.46275 0.11373;
            0.40000 0.40000 0.40000];

        % ===== 连续/发散色带（返回抽样K色；若未给K则返回函数句柄）=====
    case {'viridis','plasma','magma','inferno'}
        C = colormap_sampled(key, K);
        return
    case 'icefire'  % 发散：蓝—白—红（冷暖均衡）
        C = diverging_icefire(K);
        return
    case 'bluered'  % 发散：深蓝—白—深红
        C = diverging_bluered(K);
        return

        % ===== 大规模分类色（算法生成）=====
    case 'categories30'
        if isempty(K), K = 30; end
        C = distinct_colors(K, 'muted');
    case 'categories20'
        if isempty(K), K = 20; end
        C = distinct_colors(K, 'bright');

    otherwise
        error('Unknown palette key: %s', which);
end

% 若用户指定K并小于行数，裁剪前K；若大于行数，做扩展
if ~isempty(K)
    if K <= size(C,1)
        C = C(1:K, :);
    else
        C = extend_palette(C, K);
    end
end

% ------------ 全局调整：变暗/增亮/去饱和 ------------
if lval > 0, C = smf_lighten(C, lval);  end
if dval > 0, C = smf_darken(C, dval);   end
if sval > 0, C = smf_desaturate(C, sval); end

% 保证合法范围
C = min(max(C,0),1);
end

% ===== 简便API：从库里抽K个（不足自动扩展）=====
function C = smf_color_pick(key, K, varargin)
C0 = smf_colorlib(key, varargin{:});
if K <= size(C0,1), C = C0(1:K,:); else, C = extend_palette(C0, K); end
end

% ====== 工具：整体加深/增亮/去饱和 ======
function C2 = smf_darken(C, d)
C2 = (1-d).*C;
tooLight = mean(C,2) > 0.94;                        % 极浅保护
C2(tooLight,:) = max(0.90, 1-d).*C(tooLight,:);
end
function C2 = smf_lighten(C, l)
C2 = C + (1-C).*l;
end
function C2 = smf_desaturate(C, s)
% 线性插值到灰（相同亮度近似）：C*(1-s) + gray*s
g = mean(C,2); G = [g g g];
C2 = C*(1-s) + G*s;
end

% ====== 扩展调色板到K：基色循环 + 亮度/饱和度微调 ======
function Ck = extend_palette(C, K)
n = size(C,1);
Ck = zeros(K,3);
for i=1:K
    base = C( mod(i-1,n)+1 , : );
    cyc  = floor((i-1)/n);             % 第几轮
    % 每轮交替做轻微的亮度/饱和度扰动（不破坏主色相）
    delta = 0.10 * ((-1)^(cyc)) * min(1, 0.6 + 0.1*cyc);
    if mod(cyc,2)==0
        tmp = smf_lighten(base, max(0, 0.08*cyc));
    else
        tmp = smf_darken(base,  max(0, 0.06*cyc));
    end
    tmp = smf_desaturate(tmp, max(0, 0.05*cyc));
    Ck(i,:) = min(max(tmp + delta*[0 0 0],0),1);
end
end

% ====== 大规模分类色生成器（HSV等距 + 亮度交错；muted更稳）=====
function C = distinct_colors(K, mode)
if nargin<2, mode='muted'; end
phi = (1+sqrt(5))/2;                 % 黄金比
h = mod((0:K-1)*(1/phi), 1);         % 均匀铺色相
switch lower(mode)
    case 'bright', S = 0.70; V1 = 0.90; V2 = 0.75;
    otherwise,    S = 0.55; V1 = 0.85; V2 = 0.70; % muted：更适合大量类别
end
C = zeros(K,3);
for i=1:K
    v = V1; if mod(i,2)==0, v = V2; end
    C(i,:) = hsv2rgb([h(i), S, v]);
end
% 避免过浅：整体小幅加深
C = smf_darken(C, 0.05);
end

% ====== 连续色带抽样（viridis/plasma/magma/inferno）=====
function C = colormap_sampled(name, K)
if isempty(K)
    C = @(k) colormap_sampled(name, k); % 返回函数句柄，兼容用法
    return
end
t = linspace(0,1,K).';
switch name
    case 'viridis', C = viridis_core(t);
    case 'plasma',  C = plasma_core(t);
    case 'magma',   C = magma_core(t);
    case 'inferno', C = inferno_core(t);
end
end

% --- 近似核（基于公开插值系数的光滑近似；0-1输入t） ---
function C = viridis_core(t)
C = [ ...
    0.267 0.004 0.329;
    0.283 0.141 0.458;
    0.254 0.265 0.530;
    0.207 0.372 0.553;
    0.164 0.471 0.558;
    0.128 0.567 0.551;
    0.135 0.659 0.518;
    0.266 0.748 0.441;
    0.477 0.821 0.318;
    0.741 0.873 0.150;
    0.993 0.906 0.144];
C = interp1(linspace(0,1,size(C,1))', C, t, 'pchip');
end
function C = plasma_core(t)
C = [ ...
    0.050 0.030 0.527;
    0.283 0.008 0.641;
    0.477 0.030 0.664;
    0.647 0.116 0.620;
    0.791 0.238 0.554;
    0.902 0.384 0.488;
    0.973 0.554 0.425;
    0.996 0.736 0.369;
    0.967 0.906 0.328];
C = interp1(linspace(0,1,size(C,1))', C, t, 'pchip');
end
function C = magma_core(t)
C = [ ...
    0.001 0.000 0.014;
    0.063 0.043 0.231;
    0.200 0.120 0.429;
    0.360 0.176 0.518;
    0.525 0.235 0.513;
    0.682 0.306 0.471;
    0.820 0.407 0.397;
    0.929 0.553 0.300;
    0.992 0.743 0.235;
    0.987 0.991 0.749];
C = interp1(linspace(0,1,size(C,1))', C, t, 'pchip');
end
function C = inferno_core(t)
C = [ ...
    0.002 0.000 0.014;
    0.122 0.047 0.282;
    0.254 0.109 0.382;
    0.402 0.164 0.414;
    0.561 0.225 0.404;
    0.716 0.312 0.353;
    0.852 0.438 0.268;
    0.948 0.614 0.160;
    0.987 0.815 0.090;
    0.988 0.998 0.644];
C = interp1(linspace(0,1,size(C,1))', C, t, 'pchip');
end

% ====== 两个发散色带 ======
function C = diverging_bluered(K)
if isempty(K), K=11; end
t = linspace(0,1,K).';
blue = [0.019 0.188 0.380];
white= [1 1 1];
red  = [0.647 0.058 0.082];
C = zeros(K,3);
for i=1:K
    if t(i)<=0.5
        C(i,:) = blue + (white-blue)*(t(i)/0.5);
    else
        C(i,:) = white + (red-white)*((t(i)-0.5)/0.5);
    end
end
end
function C = diverging_icefire(K)
if isempty(K), K=11; end
t = linspace(0,1,K).';
cold=[0.059 0.365 0.647]; mid=[0.985 0.985 0.985]; warm=[0.800 0.196 0.012];
C = zeros(K,3);
for i=1:K
    if t(i)<=0.5
        C(i,:) = cold + (mid-cold)*(t(i)/0.5);
    else
        C(i,:) = mid + (warm-mid)*((t(i)-0.5)/0.5);
    end
end
end

function Ck = expand_palette_strong(B, K)
% 扩展基色B到K色：按基色循环，并在每一“轮”做明显的亮度交替与轻微饱和/色相扰动
% - 亮度交替：V × [1.00, 0.78, 0.92, 0.66] 循环，层次清楚、印刷友好
% - 饱和微调：S × [1.00, 0.90, 0.85, 0.95] 循环，避免高亮撞色
% - 色相抖动：每轮 +/− ~0.03 的 hue 偏移，打散近邻相似度

n = size(B,1);
if K <= n, Ck = B(1:K,:); return; end

% 基色转HSV
HB = rgb2hsv(B);

% 每轮的系数（可按需调整强度）
v_cycle = [1.00, 0.78, 0.92, 0.66];   % 明→暗→中→更暗（明显交替）
s_cycle = [1.00, 0.90, 0.85, 0.95];   % 适度去饱和
h_jitt  = [ +0.00, +0.03, -0.02, +0.02 ]; % 轻微色相扰动（避免撞色）

Ck = zeros(K,3);
for i = 1:K
    idx = mod(i-1, n) + 1;           % 基色索引
    cyc = floor((i-1) / n) + 1;      % 第几轮（从1开始）

    h = HB(idx,1);
    s = HB(idx,2);
    v = HB(idx,3);

    % 当前轮的系数（循环取）
    j  = mod(cyc-1, numel(v_cycle)) + 1;
    vv = v * v_cycle(j);
    ss = s * s_cycle(j);
    hh = h + h_jitt(j);

    % 轻微黄金比抖动（跨很多色时进一步打散）
    hh = hh + 0.012 * mod(cyc, 5);   % 0,0.012,0.024,0.036,0.048...
    hh = hh - floor(hh);             % wrap 到 [0,1)

    % 裁剪合法范围
    ss = min(max(ss, 0), 1);
    vv = min(max(vv, 0), 1);

    Ck(i,:) = hsv2rgb([hh, ss, vv]);
end

% 过亮保护与轻度统一加深，避免“发虚”
tooLight = mean(Ck,2) > 0.965;
Ck(tooLight,:) = (1-0.06).*Ck(tooLight,:);  % 轻压亮度
Ck = min(max(Ck, 0), 1);
end

