clc; clear; close all;

%% =========================================================
%  单腿动力学仿真与验证（对应论文 3.3 节）
%  功能：
%  1) 计算单腿一个运动周期内髋/膝关节总力矩
%  2) 分解动力学项与支撑载荷项
%  3) 绘制论文所需三张图
%
%  机器人参数：
%  机身质量 M_body  = 9.0 kg
%  单腿质量 M_leg   = 2.8 kg（不含轮）
%  单轮质量 M_wheel = 1.2 kg（含轮毂电机）
%  总质量 M_total   = 25.0 kg
%% =========================================================

%% 1. 基本参数
g = 9.81;             % 重力加速度 m/s^2

% 连杆长度
l1 = 0.30;            % 大腿长度 m
l2 = 0.30;            % 小腿长度 m

% 质量参数
M_body  = 9.0;        % 机身质量 kg
M_leg   = 2.8;        % 单腿质量 kg（不含轮）
M_wheel = 1.2;        % 单轮质量 kg（含轮毂电机）
M_total = M_body + 4*(M_leg + M_wheel);   % 整机总质量 kg = 25

% 单腿质量等效分配（与2.6.3一致：按长度比例分配）
m1 = M_leg * l1 / (l1 + l2);   % 大腿等效质量
m2 = M_leg * l2 / (l1 + l2);   % 小腿等效质量
mw = M_wheel;                  % 足端集中质量

% 转动惯量（均匀杆近似）
I1 = (1/12) * m1 * l1^2;
I2 = (1/12) * m2 * l2^2;

%% 2. 支撑载荷设定
% 对应论文3.3的验证工况：采用较保守的单腿支撑等效载荷
% 这里按三足支撑 + 安全系数1.2 进行校核
support_legs = 3;
safety_factor = 1.2;
Fd = safety_factor * (M_total * g / support_legs);   % 单腿等效竖向载荷 N

%% 3. 仿真时间
T = 2.0;               % 周期 s
N = 600;               % 采样点数
t = linspace(0, T, N);
w = 2*pi/T;            % 角频率

%% 4. 关节轨迹设置
% 用周期函数模拟支撑相与摆动相中的关节变化
% 与前文说明一致，取一个较典型的中间姿态变化范围
q2_mid_deg = 35;       % 髋屈伸中心角 deg
q2_amp_deg = 15;       % 髋屈伸振幅 deg

q3_mid_deg = 80;       % 膝屈伸中心角 deg
q3_amp_deg = 20;       % 膝屈伸振幅 deg

phi3 = pi/6;           % 膝关节相位差

% 关节角
q2 = deg2rad(q2_mid_deg + q2_amp_deg * sin(w*t));
q3 = deg2rad(q3_mid_deg + q3_amp_deg * sin(w*t + phi3));

% 角速度
dq2 = deg2rad(q2_amp_deg) * w * cos(w*t);
dq3 = deg2rad(q3_amp_deg) * w * cos(w*t + phi3);

% 角加速度
ddq2 = -deg2rad(q2_amp_deg) * w^2 * sin(w*t);
ddq3 = -deg2rad(q3_amp_deg) * w^2 * sin(w*t + phi3);

%% 5. 初始化变量
tau_dyn_hip  = zeros(1, N);    % 髋关节动力学项
tau_dyn_knee = zeros(1, N);    % 膝关节动力学项

tau_load_hip  = zeros(1, N);   % 髋关节支撑载荷项
tau_load_knee = zeros(1, N);   % 膝关节支撑载荷项

tau_total_hip  = zeros(1, N);  % 髋关节总力矩
tau_total_knee = zeros(1, N);  % 膝关节总力矩

x_foot = zeros(1, N);          % 足端x坐标
z_foot = zeros(1, N);          % 足端z坐标

%% 6. 主循环
for i = 1:N
    % 当前时刻关节状态
    q2_i = q2(i);
    q3_i = q3(i);
    dq_i = [dq2(i); dq3(i)];
    ddq_i = [ddq2(i); ddq3(i)];

    c2  = cos(q2_i);
    s2  = sin(q2_i);
    c3  = cos(q3_i);
    s3  = sin(q3_i);
    c23 = cos(q2_i + q3_i);
    s23 = sin(q2_i + q3_i);

    %% 6.1 足端位置
    x_foot(i) = l1*c2 + l2*c23;
    z_foot(i) = -l1*s2 - l2*s23;

    %% 6.2 惯性矩阵 M(q)
    M11 = I1 + I2 ...
        + m1*(l1/2)^2 ...
        + m2*(l1^2 + (l2/2)^2 + 2*l1*(l2/2)*c3) ...
        + mw*(l1^2 + l2^2 + 2*l1*l2*c3);

    M12 = I2 ...
        + m2*((l2/2)^2 + l1*(l2/2)*c3) ...
        + mw*(l2^2 + l1*l2*c3);

    M21 = M12;
    M22 = I2 + m2*(l2/2)^2 + mw*l2^2;

    Mmat = [M11, M12;
            M21, M22];

    %% 6.3 科氏/离心项 C(q,dq)
    h = m2*l1*(l2/2)*s3 + mw*l1*l2*s3;

    Cmat = [ -h*dq3(i),         -h*(dq2(i)+dq3(i));
              h*dq2(i),          0 ];

    %% 6.4 重力项 G(q)
    G1 = (m1*(l1/2) + m2*l1 + mw*l1)*g*cos(q2_i) ...
       + (m2*(l2/2) + mw*l2)*g*cos(q2_i + q3_i);

    G2 = (m2*(l2/2) + mw*l2)*g*cos(q2_i + q3_i);

    Gvec = [G1; G2];

    %% 6.5 动力学项力矩
    % tau_dyn = M(q)ddq + C(q,dq)dq + G(q)
    tau_dyn = Mmat * ddq_i + Cmat * dq_i + Gvec;
    tau_dyn_hip(i)  = tau_dyn(1);
    tau_dyn_knee(i) = tau_dyn(2);

    %% 6.6 雅可比矩阵 J(q)
    J = [ -l1*sin(q2_i) - l2*sin(q2_i + q3_i),   -l2*sin(q2_i + q3_i);
          -l1*cos(q2_i) - l2*cos(q2_i + q3_i),   -l2*cos(q2_i + q3_i) ];

    %% 6.7 支撑载荷项
    % 仅考虑足端竖向支撑力
    F = [0; -Fd];
    tau_load = J' * F;

    tau_load_hip(i)  = tau_load(1);
    tau_load_knee(i) = tau_load(2);

    %% 6.8 总力矩
    tau_total = tau_dyn + tau_load;
    tau_total_hip(i)  = tau_total(1);
    tau_total_knee(i) = tau_total(2);
end

%% 7. 峰值计算
hip_peak  = max(abs(tau_total_hip));
knee_peak = max(abs(tau_total_knee));

hip_dyn_peak  = max(abs(tau_dyn_hip));
knee_dyn_peak = max(abs(tau_dyn_knee));

hip_load_peak  = max(abs(tau_load_hip));
knee_load_peak = max(abs(tau_load_knee));

%% 8. 电机峰值参数
J60_6_peak  = 19.94;    % N·m
J60_10_peak = 30.50;    % N·m

%% 9. 输出结果
fprintf('====================================================\n');
fprintf('         单腿动力学仿真结果（对应论文3.3）\n');
fprintf('====================================================\n');
fprintf('机身质量 M_body   = %.2f kg\n', M_body);
fprintf('单腿质量 M_leg    = %.2f kg\n', M_leg);
fprintf('单轮质量 M_wheel  = %.2f kg\n', M_wheel);
fprintf('整机质量 M_total  = %.2f kg\n', M_total);
fprintf('单腿等效载荷 Fd   = %.2f N\n', Fd);

fprintf('\n---------------- 力矩峰值 ----------------\n');
fprintf('髋关节动力学项峰值     = %.2f N·m\n', hip_dyn_peak);
fprintf('膝关节动力学项峰值     = %.2f N·m\n', knee_dyn_peak);
fprintf('髋关节支撑载荷项峰值   = %.2f N·m\n', hip_load_peak);
fprintf('膝关节支撑载荷项峰值   = %.2f N·m\n', knee_load_peak);
fprintf('髋关节总力矩峰值       = %.2f N·m\n', hip_peak);
fprintf('膝关节总力矩峰值       = %.2f N·m\n', knee_peak);

fprintf('\n---------------- 与电机峰值比较 ----------------\n');
fprintf('J60-6  峰值扭矩 = %.2f N·m\n', J60_6_peak);
fprintf('J60-10 峰值扭矩 = %.2f N·m\n', J60_10_peak);

if hip_peak <= J60_10_peak
    fprintf('髋关节：J60-10 可满足该工况需求。\n');
else
    fprintf('髋关节：峰值超过 J60-10，需要说明该结果属于保守校核或极限工况。\n');
end

if knee_peak <= J60_10_peak
    fprintf('膝关节：J60-10 可满足该工况需求。\n');
else
    fprintf('膝关节：峰值超过 J60-10，需要进一步优化工况或结构参数。\n');
end

%% 10. 绘图
% 说明：
% 图中不加数值注解，避免遮挡
% 只保留曲线、标题、坐标轴和图例

%% 图1：单腿关节总力矩变化曲线
figure('Name','Total Joint Torque','Color','w');
plot(t, tau_total_hip,  'LineWidth', 2); hold on;
plot(t, tau_total_knee, 'LineWidth', 2);
yline(J60_6_peak,  '--k', 'LineWidth', 1);
yline(-J60_6_peak, '--k', 'LineWidth', 1);
yline(J60_10_peak, ':k',  'LineWidth', 1.2);
yline(-J60_10_peak, ':k', 'LineWidth', 1.2);
grid on;
xlabel('t / s');
ylabel('Torque / (N·m)');
title('单腿关节力矩变化曲线');
legend('\tau_{hip,total}', '\tau_{knee,total}', ...
       'J60-6 limit', '', 'J60-10 limit', '', ...
       'Location', 'best');

%% 图2：髋关节力矩组成
figure('Name','Hip Torque Components','Color','w');
plot(t, tau_dyn_hip,   'LineWidth', 1.8); hold on;
plot(t, tau_load_hip,  '--', 'LineWidth', 1.8);
plot(t, tau_total_hip, 'LineWidth', 2.2);
yline(J60_6_peak,  '--k', 'LineWidth', 1);
yline(-J60_6_peak, '--k', 'LineWidth', 1);
yline(J60_10_peak, ':k',  'LineWidth', 1.2);
yline(-J60_10_peak, ':k', 'LineWidth', 1.2);
grid on;
xlabel('t / s');
ylabel('Torque / (N·m)');
title('髋关节力矩组成');
legend('\tau_{dyn}', '\tau_{load}', '\tau_{total}', ...
       'J60-6 limit', '', 'J60-10 limit', '', ...
       'Location', 'best');

%% 图3：膝关节力矩组成
figure('Name','Knee Torque Components','Color','w');
plot(t, tau_dyn_knee,   'LineWidth', 1.8); hold on;
plot(t, tau_load_knee,  '--', 'LineWidth', 1.8);
plot(t, tau_total_knee, 'LineWidth', 2.2);
yline(J60_6_peak,  '--k', 'LineWidth', 1);
yline(-J60_6_peak, '--k', 'LineWidth', 1);
yline(J60_10_peak, ':k',  'LineWidth', 1.2);
yline(-J60_10_peak, ':k', 'LineWidth', 1.2);
grid on;
xlabel('t / s');
ylabel('Torque / (N·m)');
title('膝关节力矩组成');
legend('\tau_{dyn}', '\tau_{load}', '\tau_{total}', ...
       'J60-6 limit', '', 'J60-10 limit', '', ...
       'Location', 'best');

%% 11. 可选：足端轨迹图（若3.3节需要可补充）
figure('Name','Foot Trajectory','Color','w');
plot(x_foot, z_foot, 'LineWidth', 2);
grid on; axis equal;
xlabel('x / m');
ylabel('z / m');
title('单腿足端轨迹');

%% 12. 结果表（便于写论文）
ResultTable = table( ...
    ["Hip Pitch"; "Knee Pitch"], ...
    [hip_dyn_peak; knee_dyn_peak], ...
    [hip_load_peak; knee_load_peak], ...
    [hip_peak; knee_peak], ...
    'VariableNames', {'Joint', 'DynTorquePeak_Nm', 'LoadTorquePeak_Nm', 'TotalTorquePeak_Nm'} ...
    );

disp(' ');
disp('================ 结果汇总表 ================');
disp(ResultTable);