% MATLAB 脚本：利用牛顿法、不动点迭代与埃特金加速法求解方程
%    x^3 - sin(x) - 12*x + 1 = 0

% 清除环境变量
clear; clc;

% 设置容忍度和最大迭代次数
tol = 1e-6;
max_iter = 100;

%% 使用牛顿法求解
disp('使用牛顿法求解：');
x0_newton = [-3.5, 0.1, 3.5];
for x0 = x0_newton
    root = newton_method(x0, tol, max_iter);
    fprintf('%.6f\n', root);
end

%% 使用 phi_1 不动点迭代法求解（适用于区间 [-4,-3] 与 [3,4]）
disp(' ');
disp('使用 phi_1 不动点迭代法求解（适用于区间 [-4,-3] 与 [3,4]）：');
x0_phi1 = [-3.5, 3.5];
for x0 = x0_phi1
    root = fixed_point_phi1(x0, tol, max_iter);
    fprintf('%.6f\n', root);
end

%% 使用 phi_2 不动点迭代法求解（适用于区间 [0,0.2]）
disp(' ');
disp('使用 phi_2 不动点迭代法求解（适用于区间 [0,0.2]）：');
root_phi2 = fixed_point_phi2(0.1, tol, max_iter);
fprintf('%.6f\n', root_phi2);

%% 使用埃特金加速法对 phi_2 迭代序列求解（适用于区间 [0,0.2]）
disp(' ');
disp('使用埃特金加速法对 phi_2 迭代序列求解（适用于区间 [0,0.2]）：');
% 调用时传入 false，不打印内部提示信息
root_aitken = aitken_acceleration(@fixed_point_phi2, 0.1, tol, max_iter);
fprintf('%.6f\n', root_aitken);

%% --- 以下为局部函数定义 ---

function y = f(x)
    % 目标函数 f(x) = x^3 - sin(x) - 12*x + 1
    y = x.^3 - sin(x) - 12*x + 1;
end

function y = df(x)
    % f(x) 的导数：df/dx = 3*x^2 - cos(x) - 12
    y = 3*x.^2 - cos(x) - 12;
end

function root = newton_method(x0, tol, max_iter)
    % 牛顿法求根，当迭代达到容忍度时打印迭代步数并返回当前近似值
    x = x0;
    for i = 1:max_iter
        fx = f(x);
        dfx = df(x);
        if abs(dfx) < 1e-12  % 避免除以接近零的数
            fprintf('牛顿法在第 %d 步遇到导数过小，提前退出\n', i);
            root = x;
            return;
        end
        x_new = x - fx/dfx;
        if abs(x_new - x) < tol
            fprintf('牛顿法在第 %d 步达到容忍度条件\n', i);
            root = x_new;
            return;
        end
        x = x_new;
    end
    root = x;
end

function y = cube_root(val)
    % 计算立方根（对负数也正确）
    y = sign(val) .* abs(val).^(1/3);
end

function root = fixed_point_phi1(x0, tol, max_iter, verbose)
    % 不动点迭代法：φ₁(x) = (12*x + sin(x) - 1)^(1/3)
    % 可选参数 verbose 控制是否打印提示信息，默认为 true
    if nargin < 4
        verbose = true;
    end
    x = x0;
    for i = 1:max_iter
        val = 12*x + sin(x) - 1;
        x_new = cube_root(val);
        if abs(x_new - x) < tol
            if verbose
                fprintf('phi_1 不动点法在第 %d 步达到容忍度条件\n', i);
            end
            root = x_new;
            return;
        end
        x = x_new;
    end
    root = x;
end

function root = fixed_point_phi2(x0, tol, max_iter, verbose)
    % 不动点迭代法：φ₂(x) = (x^3 - sin(x) + 1) / 12
    % 可选参数 verbose 控制是否打印提示信息，默认为 true
    if nargin < 4
        verbose = true;
    end
    x = x0;
    for i = 1:max_iter
        x_new = (x^3 - sin(x) + 1) / 12;
        if abs(x_new - x) < tol
            if verbose
                fprintf('phi_2 不动点法在第 %d 步达到容忍度条件\n', i);
            end
            root = x_new;
            return;
        end
        x = x_new;
    end
    root = x;
end

function root = aitken_acceleration(phi, x0, tol, max_iter)
    % 埃特金加速法，对不动点迭代序列进行加速
    % 在内部调用 phi 时，不打印提示信息（verbose = false）
    x0_ait = x0;
    for i = 1:max_iter
        x1 = phi(x0_ait, tol, max_iter, false);
        x2 = phi(x1, tol, max_iter, false);
        denominator = x2 - 2*x1 + x0_ait;
        if abs(denominator) < 1e-12
            fprintf('埃特金加速法在第 %d 步分母过小，提前退出\n', i);
            root = x2;
            return;
        end
        x_acc = x0_ait - ((x1 - x0_ait)^2) / denominator;
        if abs(x_acc - x0_ait) < tol
            fprintf('埃特金加速法在第 %d 步达到容忍度条件\n', i);
            root = x_acc;
            return;
        end
        x0_ait = x_acc;
    end
    root = x0_ait;
end
