% ***      ——       差分变换制导算法       ——      *** %
% ***  Differential Transformation Guidance  *** %
% ***       月球软着陆燃料最优和能量最优控制       *** %

clc; clear; close all;


% 仿真参数
Ts = 0.5;                % 采样时间/仿真步长 (s)
dt = 0.01;              % 数值积分步长
n_taylor = 15;       % 泰勒展开项数 n
fuel_opti = true;   % 优化问题：true燃料最优，false能量最优
sim_case = 2;         % 仿真情况：1为90°倾角，2为45°倾角


% 月球参数
global u Radius g_Isp
u = 4903 * 1e9;            % 引力常数 (m3/s2)
Radius = 1738 * 1e3;   % 月球半径 (m)


% 着陆器参数
g_Isp = 9.81 * 315;  % 喷气速度 (m/s)
Tmax = 2200;          % 最大推力 (kg m/s2)
kmin = 0.4;              % 节流阀 kmin<=k<=1

m0 = 874.4;                     % 进入质量 (kg)
pm0 = 0;                          % m的初始协态猜测（pm0<1）
M1_P = m0 * (1-pm0);    % 常数项


% 椭圆轨道参数（进入点）
o.ra = Radius + 100 * 1e3;              % 远地点 (m)
o.rp = Radius + 15 * 1e3;                % 近地点 (m)
o.O = deg2rad(0);                            % 升交点赤经Ω
o.w = deg2rad(0);                            % 近地点幅角ω
o.theta = deg2rad(0);                      % 真近点角θ
if sim_case == 1                               % 轨道倾角i
    o.i = deg2rad(90);                        
else
    o.i = deg2rad(45);
end

o.a = (o.ra + o.rp) / 2;                      % 半长轴 (m)
o.e = (o.ra - o.rp) / (o.ra + o.rp);     % 偏心率 (1)
o.T = 2 * pi * sqrt(o.a^3 / u);           % 轨道周期 (s)
o.p = o.a * (1 - o.e^2);                      % 半通径 (m)
o.h = sqrt(o.p * u);                            % 角动量 (m2/s)
o.Va = o.h / o.ra;                               % 远地点速度 (m/s)
o.Vp = o.h / o.rp;                               % 近地点速度 (m/s)

o.r = (o.h^2/u)/(1+o.e*cos(o.theta));     % 当前轨道半径 (m)
o.V = sqrt(-u/o.a+2*u/o.r);                      % 当前轨道速度 (m/s)
o.gamma = atan(o.e*sin(o.theta)/(1+o.e*cos(o.theta))); % flight path angle

[o.r_mat, o.V_mat, o.LLA, o.orbit] = Six2Ine(o);  % 当前位置速度经纬高


% 目标参数（着陆点）
f.V_mat = [0; 0; 0];                    % 速度 =0
f.V = norm(f.V_mat);

if sim_case == 1                       % 倾角 =90
    f.Lon = deg2rad(0);              % 经度 =0                
    f.Lat = deg2rad(16.1508);    % 纬度 =16

else                                            % 倾角 =45  
    f.Lon = deg2rad(11.46113);  % 经度 =11
    f.Lat = deg2rad(11.46146);   % 纬度 =11
end
f.Alt = 0;                                     % 高度 =0

f.r = Radius + f.Alt;   
f.rx = f.r * cos(f.Lat) * cos(f.Lon);
f.ry = f.r * cos(f.Lat) * sin(f.Lon);
f.rz = f.r * sin(f.Lat);
f.r_mat = [f.rx; f.ry; f.rz];


% 初始状态
t0 = 0;

r0_mat = o.r_mat;
rx0 = r0_mat(1);
ry0 = r0_mat(2);
rz0 = r0_mat(3);
r0 = o.r;

V0_mat = o.V_mat;
Vx0 = V0_mat(1);
Vy0 = V0_mat(2);
Vz0 = V0_mat(3);

Lon0 = o.LLA(1);
Lat0 = o.LLA(2);
Alt0 = o.LLA(3);

R = eye(6,6);  tGo0 = 1;  % 防止R和tGo奇异用  


% 两点边值条件
x0_bar = [ r0_mat; V0_mat ];
xf_bar = [ f.r_mat; f.V_mat ];


% 目标北东地坐标系
NEG = [ -cos(f.Lon)*sin(f.Lat),   -sin(f.Lon)*sin(f.Lat),    cos(f.Lat);
             -sin(f.Lon),                   cos(f.Lon),                  0;             
             -cos(f.Lon)*cos(f.Lat),  -sin(f.Lon)*cos(f.Lat),  -sin(f.Lat);  ];
bias = f.r_mat;


 

% *****************缓存库******************* %
data.t = [t0];       
data.rx = [rx0];    % 状态(月心惯性系)
data.ry = [ry0];
data.rz = [rz0];
data.vx = [Vx0];
data.vy = [Vy0];
data.vz = [Vz0];
data.m = [m0];
data.pm = [pm0];  % 质量协态

data.altitude = [Alt0];       % 高度
data.longitude = [Lon0];   % 经度
data.latitude = [Lat0];       % 纬度

data.tgo = [];               % 剩余时间
data.downrange = [];  % 剩余航程
data.crossrange = [];  % 剩余横程
data.xyz = [(NEG*(r0_mat-bias))'];  % 目标北东地坐标

data.k = [];     % 油门大小
data.T = [];     % 推力大小
data.iT = [];    % 推力方向
data.u = [];     % 控制量(推力加速度大小）
data.aT = [];   % 推力方向(最大推力加速度大小）
% ********************************************* %





%%%   ——   制导系统   ——   %%%   


%  ###  DT制导算法  ### %
tic
while true
    % —— 计算当前tGo —— %
    [tGo1, D_R0, C_R0] = GetTgo(V0_mat, Lon0, Lat0, r0, f);
    if norm(V0_mat) <= 0
        tGo1 = tGo0;    % 奇异用上一时刻的值
    end
    tGo0 = tGo1;
    
    % —— 计算当前协态 —— %
    % 计算A矩阵
    A1 = [  zeros(3,3),                            eye(3,3);
                -u/r0^3 * eye(3,3),           zeros(3,3);  ];

    A2_2 = u/r0^5 * [  (ry0^2+rz0^2-2*rx0^2),    (-3*rx0*ry0),    (-3*rx0*rz0);
                                    (-3*rx0*ry0),    (rx0^2+rz0^2-2*ry0^2),    (-3*rz0*ry0);
                                    (-3*rx0*rz0),    (-3*rz0*ry0),    (ry0^2+rx0^2-2*rz0^2);  ];
    
    A2 = [  zeros(3,3),  A2_2;
                -eye(3,3),    zeros(3,3);  ];

    B1 = [  zeros(3,3),  zeros(3,3);
                zeros(3,3),  -eye(3,3);    ];
    if fuel_opti   %燃料最优
        %B1 = g_Isp * Tmax/ (m0^2*(1-pm0)) * B1;
        B1 = g_Isp * Tmax/ (m0 * M1_P) * B1;
    end
    
    A = [  A1,              B1;
              zeros(6,6),  A2;  ];
    
    % 计算A的j次方，迭代计算Q和R
    Aj = eye(12,12);        % A的0次方
    Qj = x0_bar;              % Q的0级数
    Rj = zeros(6,6);         % R的0级数

    for j =1:1:n_taylor
        Aj = Aj * A;               % 从A^1到A^n
        Aj1 = Aj(1:6,1:6);
        Aj2 = Aj(1:6,7:12);

        Qj = Qj + tGo0^j / factorial(j) * Aj1 * x0_bar;
        Rj = Rj + tGo0^j / factorial(j) * Aj2;
    end
    %if det(Rj)>0
    %R = Rj;  % R奇异时用上一步的
    %end
    
    % 协态变量
    %p0 = R \ ( xf_bar - Qj );
    p0 = Rj \ ( xf_bar - Qj );
   
    % —— 计算当前控制量 —— %
    pv = p0(4:6,1);

    if fuel_opti  % 燃料最优
        %S = (m0*(1-pm0))/g_Isp - g_Isp/(m0*(1-pm0)) * norm(pv)^2;
        S = M1_P/g_Isp - g_Isp/M1_P * norm(pv)^2;
        if S>0
            k = kmin;
        else
            k = 1;
        end
        %aT = - g_Isp * Tmax / (m0^2*(1-pm0)) * pv;   
        aT = - g_Isp * Tmax / (m0*M1_P) * pv;   
        control = k * aT;
        T = k * Tmax;
        dpm = norm(control) * pm0 / g_Isp;
    else             % 能量最优
        control = -pv;
        T = norm(control)*m0;
        if T > Tmax || T < kmin*Tmax
             control = Tmax / (m0*norm(pv)) * control;        
             T = norm(control)*m0;
        end
        aT = control * Tmax / T;
        dpm = norm(control) * pm0 * Tmax/ g_Isp * (T+1e-8);
    end

    control_bar = [zeros(3,1); control];

    [iT.ox,iT.xoy,~] = Ine2LLA(aT(1),aT(2),aT(3)); %推力方向

    % ********************************************* %
    data.k = [data.k; T/Tmax];
    data.T = [data.T; T];
    data.iT = [data.iT; [iT.ox,iT.xoy] ];
    data.u = [data.u; control' ];
    data.aT = [data.aT; aT' ];
    data.tgo = [data.tgo; tGo0];               
    data.downrange = [data.downrange; D_R0]; 
    data.crossrange = [data.crossrange; C_R0]; 
    % ********************************************* %


    % —— 是否结束制导 —— %
    if r0 - f.r <= 1 || m0 <= 10
        break
    end

    % —— 获取新状态 —— %
    for ti = 0 : dt : Ts
        %[x0_bar, m0] = Runge_Kutta(x0_bar, m0, control_bar, dt);   % 四阶龙格库塔积分
        [x0_bar, m0] = Euler_PCM(x0_bar, m0, control_bar, dt);       % 欧拉预估校正积分
        
        pm0 = pm0 + dt * dpm;
    end
    
    rx0 = x0_bar(1);
    ry0 = x0_bar(2);
    rz0 = x0_bar(3);
    r0 = sqrt(rx0^2 + ry0^2 + rz0^2);

    Vx0 = x0_bar(4);
    Vy0 = x0_bar(5);
    Vz0 = x0_bar(6);
    V0_mat = [Vx0;Vy0;Vz0];
    
    % —— 获取新经纬高 —— %
    [Lon0, Lat0, Alt0] = Ine2LLA(rx0, ry0, rz0);

    % —— 获取新时间 —— %
    t0 = t0 + Ts;


    % ********************************************* %
    data.t = [data.t;t0];       
    data.rx = [data.rx;rx0];    
    data.ry = [data.ry;ry0];
    data.rz = [data.rz;rz0];
    data.vx = [data.vx;Vx0];
    data.vy = [data.vy;Vy0];
    data.vz = [data.vz;Vz0];
    data.m = [data.m;m0];
    data.pm = [data.pm;pm0];
    data.altitude = [data.altitude;Alt0];       
    data.longitude = [data.longitude;Lon0];   
    data.latitude = [data.latitude;Lat0];     
    data.xyz = [data.xyz; (NEG*(x0_bar(1:3)-bias))' ];
    % ********************************************* %
end
toc




% ###  估计 Tgo 、航程 、横程   ### %
function [tgo, D_R, C_R] = GetTgo(V_mat, Lon, Lat, r, f)
    global Radius
    % 地面系速度
    Vg_mat = Ine2Neg(V_mat, Lon, Lat);
    
    % 方位角
    Az = atan2(Vg_mat(2), Vg_mat(1)); 

    % 航向角
    c1 = cos(f.Lat) * sin(Lon - f.Lon);
    c2 = cos(Lat) * sin(f.Lat);
    c3 = sin(Lat) * cos(f.Lat) * cos(Lon - f.Lon);
    Cmat = [c1; c2; c3];

    %sRhoC = -cos(Az)*c1 + sin(Az)*c2 - sin(Az)*c3;
    sRhoC = [-cos(Az), sin(Az), -sin(Az)] * Cmat;
    RhoC = asin(sRhoC);

    %sRhoD = ( sin(Az)*c1 + cos(Az)*c2 - cos(Az)*c3 ) /  ( cos(RhoC) + 1e-10 );
    sRhoD = [sin(Az), cos(Az), -cos(Az)] * Cmat / cos(RhoC); %可能奇异
    RhoD = asin(sRhoD);
    
    % 平均距离
    Rave = ( r + f.r ) / 2;

    % 航向和横向距离
    D_R = Rave * RhoD;
    C_R = Rave * RhoC;

    % 待飞时间
    tgo = 2 * sqrt( (r-Radius)^2 + D_R^2 + C_R^2 ) / ( norm(V_mat) + f.V ); %可能奇异
end





%%%   ——   导航系统   ——   %%%


% 轨道六要素 -> 惯性系位置速度
function [r, V, LLA, orbit] = Six2Ine(o)
    global u 
    r_pqw = o.r * [  cos(o.theta);
                               sin(o.theta);
                               0;                     ];
    V_pqw = u/o.h * [  -sin(o.theta);
                                    o.e + cos(o.theta);
                                    0;                              ];

    Q1 = [  cos(o.w),   -sin(o.w),   0;
                 sin(o.w),    cos(o.w),   0;
                 0,                0,               1;  ];
    Q2 = [  1,  0,              0;
                 0,  cos(o.i),   -sin(o.i);
                 0,  sin(o.i),    cos(o.i);  ];
    Q3 = [  cos(o.O),  -sin(o.O),  0;
                 sin(o.O),   cos(o.O),  0;
                 0,              0,              1;];
    Q = Q3 * Q2 * Q1;
    
    r = Q * r_pqw;
    V = Q * V_pqw;

    [Lon, Lat, Alt] = Ine2LLA(r(1), r(2), r(3));
    LLA = [Lon; Lat; Alt];
    
    theta = -pi:pi/100:pi;  % 1*n
    ri = o.h^2/u./(1+o.e*cos(theta)); % 1*n
    orbit_pqw = [ri.*cos(theta); ri.*sin(theta); ri.*0]; % 3*n
    orbit = Q * orbit_pqw; % 3*n
    orbit = orbit'; % n*3

end

% 惯性系 -> 北东地坐标系
function Vg = Ine2Neg(V, Lon, Lat)
    T = [ -cos(Lon)*sin(Lat),   -sin(Lon)*sin(Lat),    cos(Lat);
             -sin(Lon),                   cos(Lon),                  0;             
             -cos(Lon)*cos(Lat),  -sin(Lon)*cos(Lat),  -sin(Lat);  ];
    Vg = T * V;  
end

% 惯性系 -> 经纬高
function [Lon, Lat, Alt] = Ine2LLA(x, y, z)
    global Radius
    Lon = atan2(y, x);
    Lat = atan( z / sqrt( x^2 + y^2 ) );
    Alt = sqrt(x^2 + y^2 + z^2) - Radius;
end





%%%   ——   动力学系统   ——   %%%


% 状态转移矩阵
function A1 = GetA1(r0)
    global u
    A1 = [  zeros(3,3),                            eye(3,3);
                -u/r0^3 * eye(3,3),           zeros(3,3);  ];
end

% 欧拉预估校正积分
function [X1, m1] = Euler_PCM(X0, m0, U, dt)
    global g_Isp
    r0 = norm(X0(1:3,1));
    A0 = GetA1(r0);

    X10 = X0 + dt * ( A0 * X0 + U );

    r10 = norm(X10(1:3,1));
    A10 = GetA1(r10);

    X1 = X0 + dt/2 * ( A0 * X0 + U + A10 * X10 + U );

    m10 = m0 + dt * ( -m0 * norm(U) / g_Isp );
    m1 = m0 + dt/2 * ( -m0 * norm(U) / g_Isp  -m10 * norm(U) / g_Isp );

end

% 四级四阶龙格库塔积分
function [X1, m1] = Runge_Kutta(X0, m0, U, dt)
    global g_Isp
    X = X0;
    r = norm(X(1:3,1));
    A = GetA1(r);
    K1 = A * X + U;

    X = X + dt/2 * K1;
    r = norm(X(1:3,1));
    A = GetA1(r);
    K2 = A * X + U;

    X = X + dt/2 * K2;
    r = norm(X(1:3,1));
    A = GetA1(r);
    K3 = A * X + U;

    X = X + dt * K3;
    r = norm(X(1:3,1));
    A = GetA1(r);
    K4 = A * X + U;

    X1 = X0 + dt/6 * (K1 + 2*K2 + 2*K3 + K4);
    
    m = m0;
    K1 = -m * norm(U) / g_Isp;

    m = m0 + dt/2 * K1;
    K2 = -m * norm(U) / g_Isp;

    m = m0 + dt/2 * K2;
    K3 = -m * norm(U) / g_Isp;

    m = m0 + dt * K3;
    K4 = -m * norm(U) / g_Isp;

    m1 = m0 + dt/6 * (K1 + 2*K2 + 2*K3 + K4);

end