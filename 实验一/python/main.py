import pandas as pd
import math

def method1(a, n):
    if n == 0:
        return math.log((a + 1) / a)
    else:
        return -a * method1(a, n - 1) + 1 / n

def method2(a, n):
    if a >= n / (n + 1):
        return (2 * a + 1) / (2 * a * (a + 1) * (n + 1))
    else:
        return 0.5 * (1 / ((a + 1) * (n + 1)) + 1 / n)

# 创建数据框架，仅包含I_n列
data = {'I_n': list(range(1, 11))}

a = 0.05
df_min = pd.DataFrame(data)
# 使用apply方法为每个索引计算method1和method2的值
df_min['method_1'] = df_min['I_n'].apply(lambda n: method1(a, n)).astype(float)
df_min['method_2'] = df_min['I_n'].apply(lambda n: method2(a, n)).astype(float)

a = 15
df_max = pd.DataFrame(data)
df_max['method_1'] = df_max['I_n'].apply(lambda n: method1(a, n)).astype(float)
df_max['method_2'] = df_max['I_n'].apply(lambda n: method2(a, n)).astype(float)

df_min.to_csv("5e-2.csv", encoding='utf-8', index=False)
df_max.to_csv("15.csv", encoding='utf-8', index=False)

print(df_max)
print('------------------')
print(df_min)
