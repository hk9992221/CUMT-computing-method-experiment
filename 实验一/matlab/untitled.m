% 清空环境变量
clear; clc;

% 定义 I_n 列（1~10），转换为列向量
I_n = (1:10)';

%% 针对 a = 0.05 的计算
a = 0.05;
method_1_min = zeros(10,1);
method_2_min = zeros(10,1);
for i = 1:length(I_n)
    n = I_n(i);
    method_1_min(i) = method1(a, n);
    method_2_min(i) = method2(a, n);
end

% 将结果存放在表格中（各列均为 double 类型）
df_min = table(I_n, method_1_min, method_2_min, ...
    'VariableNames', {'I_n','method_1','method_2'});

%% 针对 a = 15 的计算
a = 15;
method_1_max = zeros(10,1);
method_2_max = zeros(10,1);
for i = 1:length(I_n)
    n = I_n(i);
    method_1_max(i) = method1(a, n);
    method_2_max(i) = method2(a, n);
end

% 将结果存放在表格中（各列均为 double 类型）
df_max = table(I_n, method_1_max, method_2_max, ...
    'VariableNames', {'I_n','method_1','method_2'});

%% 显示结果
disp('a = 0.05')
disp(df_min)
disp('--------------')
disp('a = 15')
disp(df_max)

%% 提取变量并保存到 exp1.mat 文件
% 变量 A 为 df_min 的 method_1，B 为 df_min 的 method_2
% 变量 C 为 df_max 的 method_1，D 为 df_max 的 method_2
A = double(df_min.method_1);
B = double(df_min.method_2);
C = double(df_max.method_1);
D = double(df_max.method_2);

save('exp1.mat', 'A', 'B', 'C', 'D');

%%写入csv文件
writetable(df_min, '5e-2.csv');
writetable(df_max, '15.csv');
%% 辅助函数
function out = method1(a, n)
    if n == 0
        out = log((a + 1) / a);
    else
        out = -a * method1(a, n - 1) + 1 / n;
    end
end

function out = method2(a, n)
    if a >= n / (n + 1)
        out = (2 * a + 1) / (2 * a * (a + 1) * (n + 1));
    else
        out = 0.5 * (1 / ((a + 1) * (n + 1)) + 1 / n);
    end
end
