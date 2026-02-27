clc; clear; close all;

%% =========================================================
%  杆件长度设计：扫描 r=L2/L1，以最大三维体积 V_alpha 选最优
%  坐标系约定：X向前，Y向左，Z向上；腿自然下垂方向为 -Z
%  关节定义：
%    q1：髋外展/内收，绕 X 轴旋转（Hip Abduction about X）
%    q2：髋俯仰，绕 Y 轴旋转（Hip Pitch about Y）
%    q3：膝俯仰，绕 Y 轴旋转（Knee Pitch about Y）
%% =========================================================
% 由步幅/抬腿高度反推总腿长
S  = 0.25;      % 最大水平可达距离/步幅相关 (m)
Hf = 0.15;      % 最大抬脚高度 (m)
k  = 1.35;      % 安全系数

Lmin   = sqrt(S^2 + Hf^2);
Ltotal = k * Lmin;     % 总腿长 L1+L2

% 扫描比例范围：r = L2/L1
r_list = linspace(0.6, 1.6, 41);

% 关节范围（度）
q1_lim_deg = [-35 35];     % q1：绕X
q2_lim_deg = [-110 70];    % q2：绕Y
q3_lim_deg = [10 140];     % q3：绕Y（避开奇异/极限）

% 扫描采样密度
n1s = 13; n2s = 21; n3s = 21;

% 最优解可视化采样密度
n1 = 25;  n2 = 35;  n3 = 35;

% Alpha Shape 参数
alpha = 0.060;

% 动画参数（闭合轨迹扫边界）
T_total  = 6;     % 一个闭合周期（秒）
fps      = 60;    % 帧率
usePhase = true;  % 加相位差

%% ===============【扫描 r 计算 V_alpha(r)】==============
q1s = linspace(deg2rad(q1_lim_deg(1)), deg2rad(q1_lim_deg(2)), n1s);
q2s = linspace(deg2rad(q2_lim_deg(1)), deg2rad(q2_lim_deg(2)), n2s);
q3s = linspace(deg2rad(q3_lim_deg(1)), deg2rad(q3_lim_deg(2)), n3s);

Va = nan(size(r_list));

for ii = 1:numel(r_list)
    r  = r_list(ii);
    L1 = Ltotal/(1+r);
    L2 = r*L1;

    % 点云（解析式，q1绕X）
    % 平面（q1=0）：
    %   x  = L1*cos(q2) + L2*cos(q2+q3)
    %   z0 = -(L1*sin(q2) + L2*sin(q2+q3))
    % 绕X旋转：
    %   y = sin(q1)*z0
    %   z = cos(q1)*z0
    N = numel(q1s)*numel(q2s)*numel(q3s);
    pts = zeros(N, 3);
    idx = 0;

    for i = 1:numel(q1s)
        s1 = sin(q1s(i));
        c1 = cos(q1s(i));
        for j = 1:numel(q2s)
            q2 = q2s(j);
            cq2 = cos(q2); sq2 = sin(q2);
            for k3 = 1:numel(q3s)
                q3 = q3s(k3);
                c23 = cos(q2+q3); s23 = sin(q2+q3);

                idx = idx + 1;
                x  = L1*cq2 + L2*c23;
                z0 = -(L1*sq2 + L2*s23);
                y  = s1 * z0;
                z  = c1 * z0;
                pts(idx,:) = [x y z];
            end
        end
    end

    try
        shp = alphaShape(pts, alpha);
        Va(ii) = volume(shp);
    catch
        Va(ii) = NaN;
    end
end

% 选最优：V_alpha 最大
[Va_best, idBest] = max(Va);
r_best  = r_list(idBest);
L1_best = Ltotal/(1+r_best);
L2_best = r_best*L1_best;

%% ===============【命令行输出】==============
fprintf('\n=========== 杆件长度反推===========\n');
fprintf('输入：S=%.3f m, Hf=%.3f m, k=%.2f\n', S, Hf, k);
fprintf('Lmin=%.4f m, Ltotal=%.4f m\n', Lmin, Ltotal);
fprintf('最优（V_alpha 最大）：r_best=%.3f\n', r_best);
fprintf('L1=%.4f m, L2=%.4f m, L1+L2=%.4f m\n', L1_best, L2_best, Ltotal);
fprintf('V_alpha_max=%.6e m^3, alpha=%.3f\n', Va_best, alpha);
fprintf('===================================\n\n');

%% ===============【图1：扫描曲线】==============
figure('Name','图1 扫描曲线','Position',[120 80 900 620]);
plot(r_list, Va, 'LineWidth',2);
grid on;
xlabel('r = L_2/L_1');
ylabel('V_{alpha} (m^3)');
title('图1：三维工作空间体积 V_{alpha} 随 r 的变化');

%% ===============【表1：扫描结果表】==============
Tscan = table(r_list(:), Va(:), ...
    'VariableNames', {'r_L2除以L1','V_alpha_m3'});

figure('Name','表1 扫描结果表','Position',[1050 80 520 620]);
uitable('Data', Tscan{:,:}, ...
        'ColumnName', Tscan.Properties.VariableNames, ...
        'Units','normalized','Position',[0 0 1 1]);

%% ===============【最优 r 下的 alpha 边界】==============
q1 = linspace(deg2rad(q1_lim_deg(1)), deg2rad(q1_lim_deg(2)), n1);
q2 = linspace(deg2rad(q2_lim_deg(1)), deg2rad(q2_lim_deg(2)), n2);
q3 = linspace(deg2rad(q3_lim_deg(1)), deg2rad(q3_lim_deg(2)), n3);

N = numel(q1)*numel(q2)*numel(q3);
pts = zeros(N, 3);
idx = 0;

for i = 1:numel(q1)
    s1 = sin(q1(i));
    c1 = cos(q1(i));
    for j = 1:numel(q2)
        q2j = q2(j);
        cq2 = cos(q2j); sq2 = sin(q2j);
        for k3 = 1:numel(q3)
            q3k = q3(k3);
            c23 = cos(q2j+q3k); s23 = sin(q2j+q3k);

            idx = idx + 1;
            x  = L1_best*cq2 + L2_best*c23;
            z0 = -(L1_best*sq2 + L2_best*s23);
            y  = s1 * z0;
            z  = c1 * z0;
            pts(idx,:) = [x y z];
        end
    end
end

shp_best = alphaShape(pts, alpha);
V_alpha_best = volume(shp_best);

figure('Name','图2 Alpha边界','Position',[150 140 900 620]);
plot(shp_best, 'FaceAlpha',0.15, 'EdgeAlpha',0.30);
axis equal; grid on; view(135,25);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('图2：最优 r 下的 Alpha Shape 边界（闭合曲面）');

%% ===============【图3：坐标系与关节轴方向示意（q1绕X）】==============
figure('Name','图3 坐标系与关节轴','Position',[190 180 900 620]);
hold on; grid on; axis equal; view(135,25);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('图3：坐标系与关节轴方向示意（q1绕X，q2/q3绕Y）');

s = 0.18;  % 箭头长度
quiver3(0,0,0, s,0,0,'LineWidth',2);  % X
quiver3(0,0,0, 0,s,0,'LineWidth',2);  % Y
quiver3(0,0,0, 0,0,s,'LineWidth',2);  % Z
quiver3(0,0,0, s*1.2,0,0,'LineWidth',3);        % q1轴 ~ X
quiver3(0,0,0, 0,s*1.2,0,'LineWidth',3);        % q2轴 ~ Y
quiver3(0,0,-0.05, 0,s*1.2,0,'LineWidth',3);    % q3轴 ~ Y
axis([-0.25 0.35 -0.25 0.35 -0.35 0.25]);

%% ===============【图4：模型运动】==============
hasRTB = (exist('SerialLink','class')==8) && (exist('Link','class')==8);
leg = [];
if hasRTB
    try
        Lnk(1) = Link('revolute','d',0,'a',0 ,'alpha', 0);
        Lnk(2) = Link('revolute','d',0,'a',L1_best,'alpha', 0);
        Lnk(3) = Link('revolute','d',0,'a',L2_best,'alpha', 0);
        leg = SerialLink(Lnk, 'name', 'SingleLeg');

        leg.qlim = [deg2rad(q1_lim_deg(1)) deg2rad(q1_lim_deg(2));
                    deg2rad(q2_lim_deg(1)) deg2rad(q2_lim_deg(2));
                    deg2rad(q3_lim_deg(1)) deg2rad(q3_lim_deg(2))];
    catch
        hasRTB = false;
        leg = [];
        warning('RTB失败');
    end
end

dt = 1/fps;
q1_mid = mean(deg2rad(q1_lim_deg));  q1_amp = diff(deg2rad(q1_lim_deg))/2;
q2_mid = mean(deg2rad(q2_lim_deg));  q2_amp = diff(deg2rad(q2_lim_deg))/2;
q3_mid = mean(deg2rad(q3_lim_deg));  q3_amp = diff(deg2rad(q3_lim_deg))/2;

phi2 = 0; phi3 = 0;
if usePhase
    phi2 = pi/3;
    phi3 = -pi/5;
end

figure('Name','图4 动画扫边界','Position',[230 220 900 620]);
hold on; grid on; axis equal; view(135,25);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('图4：连续闭合运动扫掠末端轨迹');

axis_lim = 1.1*Ltotal;
axis([-axis_lim axis_lim -axis_lim axis_lim -axis_lim axis_lim]);

h_traj = plot3(NaN,NaN,NaN,'LineWidth',2);
traj = zeros(floor(T_total/dt)+2, 3);
kk = 0;

if hasRTB && ~isempty(leg)
    leg.plot([q1_mid q2_mid q3_mid], 'delay',0, ...
        'workspace',[-axis_lim axis_lim -axis_lim axis_lim -axis_lim axis_lim], 'scale',0.7);
end

t = 0;
while t <= T_total + 1e-9
    q1_t = q1_mid + q1_amp*sin(2*pi*t/T_total);
    q2_t = q2_mid + q2_amp*sin(2*pi*t/T_total + phi2);
    q3_t = q3_mid + q3_amp*sin(2*pi*t/T_total + phi3);

    if hasRTB && ~isempty(leg)
        leg.animate([q1_t q2_t q3_t]);
    end

    % 末端轨迹（解析式，严格q1绕X）
    x  = L1_best*cos(q2_t) + L2_best*cos(q2_t+q3_t);
    z0 = -(L1_best*sin(q2_t) + L2_best*sin(q2_t+q3_t));
    y  = sin(q1_t)*z0;
    z  = cos(q1_t)*z0;

    kk = kk + 1;
    traj(kk,:) = [x y z];
    set(h_traj,'XData',traj(1:kk,1),'YData',traj(1:kk,2),'ZData',traj(1:kk,3));
    drawnow limitrate

    t = t + dt;
end

fprintf('闭合性检查：||p(T)-p(0)|| = %.6e\n', norm(traj(kk,:)-traj(1,:)));