%% Newton法求解 f(x)=0 及 f'(x)=0 并绘图（单图显示）
clear; clc; close all;

%% 定义函数及其导数和二阶导数
f   = @(x) 0.1*exp(x) - sin(x).^2 + 0.5;
fp  = @(x) 0.1*exp(x) - sin(2*x);      % f'(x)=0.1e^x - sin(2x)
fpp = @(x) 0.1*exp(x) - 2*cos(2*x);     % f''(x)=0.1e^x - 2cos(2x)

%% 绘制函数及其导数图像（统一在一张图上显示）
x = linspace(-10, 3, 400);
figure;
plot(x, f(x), 'b-', 'LineWidth', 1.5); hold on;
plot(x, fp(x), 'r-', 'LineWidth', 1.5);
xlabel('x'); ylabel('y');
title('函数 f(x) 与其导数 f''(x) 图像');
legend('f(x)', 'f''(x)');
grid on;

%% Newton迭代参数设置
tol = 1e-8;     % 收敛容限
maxIter = 100;  % 最大迭代次数

%% (1) 以 x0 = -5 为初值求 f(x)=0 的零点
x0_zero = -5;
root = newtonMethod(f, fp, x0_zero, tol, maxIter);
fprintf('以 x0 = %f 求 f(x)=0 的零点： x = %f\n', x0_zero, root);

%% (2) 以 x0 = -5 为初值求 f''(x)=0 对应的极值点（求 f'(x)=0）
x0_stat1 = -5;
stat1 = newtonMethodStationary(fp, fpp, x0_stat1, tol, maxIter);
fprintf('以 x0 = %f 求 f''(x)=0 的极值点： x = %f\n', x0_stat1, stat1);

%% (3) 以 x0 = 1 为初值求 f''(x)=0 对应的极值点（求 f'(x)=0）
x0_stat2 = 1;
stat2 = newtonMethodStationary(fp, fpp, x0_stat2, tol, maxIter);
fprintf('以 x0 = %f 求 f''(x)=0 的极值点： x = %f\n', x0_stat2, stat2);

%% 局部函数：Newton法求 f(x)=0
function root = newtonMethod(f, fp, x0, tol, maxIter)
    x = x0;
    for iter = 1:maxIter
        f_val = f(x);
        fp_val = fp(x);
        if abs(fp_val) < eps
            error('导数接近0，无法继续迭代。');
        end
        x_new = x - f_val / fp_val;
        if abs(x_new - x) < tol
            root = x_new;
            return;
        end
        x = x_new;
    end
    root = x;
end

%% 局部函数：Newton法求 f'(x)=0（求极值点）
function stat = newtonMethodStationary(fp, fpp, x0, tol, maxIter)
    x = x0;
    for iter = 1:maxIter
        fp_val = fp(x);
        fpp_val = fpp(x);
        if abs(fpp_val) < eps
            error('二阶导数接近0，无法继续迭代。');
        end
        x_new = x - fp_val / fpp_val;
        if abs(x_new - x) < tol
            stat = x_new;
            return;
        end
        x = x_new;
    end
    stat = x;
end
