%% LagrangeInterpolationExperiment.m
% 演示拉格朗日插值方法对两个函数进行插值，并比较等距节点与切比雪夫节点。
% 每个图窗口使用较大的画布，避免文字与图重叠；并在底部(或顶部)放置公共 legend。
% 如果需要修改第二个函数（当前定义为 f(x)=1/(1+x^2)），
% 则需要注意并修改以下几个部分：
%
% 1. 修改函数定义：
%    在文件末尾的函数 other_function 中，将表达式
%        y = 1 ./ (1 + x.^2);
%    修改为新的函数表达式。例如，如果你想改成 f(x) = exp(-x.^2)，
%    则改成:
%        y = exp(-x.^2);
%
% 2. 修改函数说明与图标题：
%    在实验二的各个子图标题（sgtitle 和子图 title）中，
%    当前显示的是 "函数 f(x)=1/(1+x^2) ..."。如果新的函数表达式不同，
%    请相应地更新这些标题，确保图示描述与修改后的函数定义一致。
%
% 3. 修改评价点的定义（x_eval2）：
%    如果新的函数定义域与原函数不同（例如不再定义在 [-5,5] 上），
%    请在实验二的 x_eval2 定义中，更新区间范围，使得评价点覆盖合适的区域。
%
% 4. 检查其他可能的引用：
%    在实验二中，所有调用 other_function 的地方（例如生成节点处计算 f_values_eq2 与 f_values_cheb2）
%    都会使用这个函数。因此在修改函数定义后，所有这些部分都将自动应用新的函数计算值，
%    但请确认其他相关注释或图例描述与函数定义保持一致。
%
% 总结：
%    主要需要修改的是：
%      - 函数 other_function 的内部表达式；
%      - 实验二的图标题、说明文字和评价点 x_eval2 的区间（如有必要）。
%
% 修改后请确保重新运行脚本，检查所有图窗口中的显示和数值计算是否如预期。

clear; clc; close all;

%% 实验一：函数 f(x)=1/(1+25*x^2) 在 [-1,1] 上
x_eval1 = linspace(-1, 1, 400);
n_values = [2, 4, 8, 16, 32];

%% (a) 等距节点插值
figure('Position',[100, 50, 1200, 600]);  % 调整画布大小
t = tiledlayout(2, ceil(length(n_values)/2), 'TileSpacing','loose', 'Padding','compact');
for i = 1:length(n_values)
    n = n_values(i);
    nodes_eq = equidistant_nodes(-1, 1, n);
    f_values_eq = runge_function(nodes_eq);
    y_interp_eq = lagrange_interp(x_eval1, nodes_eq, f_values_eq);
    
    ax = nexttile;
    plot(x_eval1, runge_function(x_eval1), 'b-', 'LineWidth', 1.5); hold on;
    plot(x_eval1, y_interp_eq, 'g-', 'LineWidth', 1.5);
    plot(nodes_eq, f_values_eq, 'ro', 'MarkerFaceColor', 'r');
    title(sprintf('n = %d (等距)', n), 'FontSize', 12);
    set(ax,'FontSize',12);  % 坐标轴字体
    hold off;
end
sgtitle('函数 f(x)=1/(1+25x^2) 在 [-1,1] 上的等距节点插值','FontSize',14,'FontWeight','bold');

% 在底部添加公共 legend
dummyAx = axes('Position',[0 0 1 1],'Visible','off'); 
hold(dummyAx, 'on');
h1 = plot(dummyAx, nan, nan, 'b-', 'LineWidth', 1.5);
h2 = plot(dummyAx, nan, nan, 'g-', 'LineWidth', 1.5);
h3 = plot(dummyAx, nan, nan, 'ro', 'MarkerFaceColor','r');
hold(dummyAx, 'off');
lg = legend(dummyAx,[h1,h2,h3], {'原函数','插值多项式','插值节点'}, ...
    'Orientation','horizontal','FontSize',12, 'Location','southoutside');
% 若对自动摆放不满意，也可改为:
% lg.Position = [0.35 0.01 0.3 0.05];  % 手动设定

%% (b) 切比雪夫节点插值
figure('Position',[100, 50, 1200, 600]);  % 调整画布大小
t = tiledlayout(2, ceil(length(n_values)/2), 'TileSpacing','loose', 'Padding','compact');
for i = 1:length(n_values)
    n = n_values(i);
    nodes_cheb = chebyshev_nodes(-1, 1, n);
    f_values_cheb = runge_function(nodes_cheb);
    y_interp_cheb = lagrange_interp(x_eval1, nodes_cheb, f_values_cheb);
    
    ax = nexttile;
    plot(x_eval1, runge_function(x_eval1), 'b-', 'LineWidth', 1.5); hold on;
    plot(x_eval1, y_interp_cheb, 'g-', 'LineWidth', 1.5);
    plot(nodes_cheb, f_values_cheb, 'ro', 'MarkerFaceColor', 'r');
    title(sprintf('n = %d (切比雪夫)', n), 'FontSize', 12);
    set(ax,'FontSize',12);
    hold off;
end
sgtitle('函数 f(x)=1/(1+25x^2) 在 [-1,1] 上的切比雪夫节点插值','FontSize',14,'FontWeight','bold');

dummyAx = axes('Position',[0 0 1 1],'Visible','off'); 
hold(dummyAx, 'on');
h1 = plot(dummyAx, nan, nan, 'b-', 'LineWidth', 1.5);
h2 = plot(dummyAx, nan, nan, 'g-', 'LineWidth', 1.5);
h3 = plot(dummyAx, nan, nan, 'ro', 'MarkerFaceColor','r');
hold(dummyAx, 'off');
lg = legend(dummyAx,[h1,h2,h3], {'原函数','插值多项式','插值节点'}, ...
    'Orientation','horizontal','FontSize',12, 'Location','southoutside');

%% 实验二：函数 f(x)=1/(1+x^2) 在 [-5,5] 上
x_eval2 = linspace(-5, 5, 400);
n_values2 = [4, 8, 16, 32];

%% (a) 等距节点插值
figure('Position',[200, 50, 1200, 600]);  % 调整画布大小
t = tiledlayout(2, ceil(length(n_values2)/2), 'TileSpacing','loose', 'Padding','compact');
for i = 1:length(n_values2)
    n = n_values2(i);
    nodes_eq2 = equidistant_nodes(-5, 5, n);
    f_values_eq2 = other_function(nodes_eq2);
    y_interp_eq2 = lagrange_interp(x_eval2, nodes_eq2, f_values_eq2);
    
    ax = nexttile;
    plot(x_eval2, other_function(x_eval2), 'b-', 'LineWidth', 1.5); hold on;
    plot(x_eval2, y_interp_eq2, 'g-', 'LineWidth', 1.5);
    plot(nodes_eq2, f_values_eq2, 'ro', 'MarkerFaceColor', 'r');
    title(sprintf('n = %d (等距)', n), 'FontSize', 12);
    set(ax,'FontSize',12);
    hold off;
end
sgtitle('函数 f(x)=1/(1+x^2) 在 [-5,5] 上的等距节点插值','FontSize',14,'FontWeight','bold');

dummyAx = axes('Position',[0 0 1 1],'Visible','off'); 
hold(dummyAx, 'on');
h1 = plot(dummyAx, nan, nan, 'b-', 'LineWidth', 1.5);
h2 = plot(dummyAx, nan, nan, 'g-', 'LineWidth', 1.5);
h3 = plot(dummyAx, nan, nan, 'ro', 'MarkerFaceColor','r');
hold(dummyAx, 'off');
lg = legend(dummyAx,[h1,h2,h3], {'原函数','插值多项式','插值节点'}, ...
    'Orientation','horizontal','FontSize',12, 'Location','southoutside');

%% (b) 切比雪夫节点插值
figure('Position',[200, 50, 1200, 600]);  % 调整画布大小
t = tiledlayout(2, ceil(length(n_values2)/2), 'TileSpacing','loose', 'Padding','compact');
for i = 1:length(n_values2)
    n = n_values2(i);
    nodes_cheb2 = chebyshev_nodes(-5, 5, n);
    f_values_cheb2 = other_function(nodes_cheb2);
    y_interp_cheb2 = lagrange_interp(x_eval2, nodes_cheb2, f_values_cheb2);
    
    ax = nexttile;
    plot(x_eval2, other_function(x_eval2), 'b-', 'LineWidth', 1.5); hold on;
    plot(x_eval2, y_interp_cheb2, 'g-', 'LineWidth', 1.5);
    plot(nodes_cheb2, f_values_cheb2, 'ro', 'MarkerFaceColor', 'r');
    title(sprintf('n = %d (切比雪夫)', n), 'FontSize', 12);
    set(ax,'FontSize',12);
    hold off;
end
sgtitle('函数 f(x)=1/(1+x^2) 在 [-5,5] 上的切比雪夫节点插值','FontSize',14,'FontWeight','bold');

dummyAx = axes('Position',[0 0 1 1],'Visible','off'); 
hold(dummyAx, 'on');
h1 = plot(dummyAx, nan, nan, 'b-', 'LineWidth', 1.5);
h2 = plot(dummyAx, nan, nan, 'g-', 'LineWidth', 1.5);
h3 = plot(dummyAx, nan, nan, 'ro', 'MarkerFaceColor','r');
hold(dummyAx, 'off');
lg = legend(dummyAx,[h1,h2,h3], {'原函数','插值多项式','插值节点'}, ...
    'Orientation','horizontal','FontSize',12, 'Location','southoutside');

%% --------- 以下为局部函数 ---------

function L = lagrange_interp_scalar(x, nodes, f_values)
    % 计算单点 x 处的拉格朗日插值值
    N = length(nodes);
    L = 0;
    for j = 1:N
        basis = 1;
        for m = 1:N
            if m ~= j
                basis = basis * (x - nodes(m)) / (nodes(j) - nodes(m));
            end
        end
        L = L + f_values(j) * basis;
    end
end

function y = lagrange_interp(x_eval, nodes, f_values)
    % 对 x_eval 数组中每个元素计算拉格朗日插值值
    y = zeros(size(x_eval));
    for k = 1:length(x_eval)
        y(k) = lagrange_interp_scalar(x_eval(k), nodes, f_values);
    end
end

function nodes = equidistant_nodes(a, b, n)
    % 在区间 [a,b] 内生成等距节点，共 n+1 个
    nodes = linspace(a, b, n+1);
end

function nodes = chebyshev_nodes(a, b, n)
    % 在区间 [a,b] 内生成切比雪夫节点，共 n+1 个
    k = 1:(n+1);
    nodes = (a+b)/2 + (b-a)/2 * cos((2*k - 1)*pi/(2*(n+1)));
end

function y = runge_function(x)
    % 函数: f(x)=1/(1+25*x^2)
    y = 1 ./ (1 + 25*x.^2);
end

function y = other_function(x)
    % 函数: f(x)=1/(1+x^2)
    y = 1 ./ (1 + x.^2);
end
