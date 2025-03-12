import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 如果环境允许，可以用 SciPy 做数值积分:
from scipy.integrate import quad


# ============ 1. 定义数值积分函数，获得参考值 ============
def numeric_I(a, n):
    """用quad对 I_n(a) = \int_0^1 x^n / (a+x) dx 做数值积分."""

    def integrand(x):
        return (x ** n) / (a + x)

    val, _ = quad(integrand, 0, 1)
    return val


# ============ 2. 题目给出的两种方法 method1 / method2 ============
def method1(a, n):
    """
    递推格式:
      I_0 = ln((a + 1)/a)
      I_n = -a * I_{n-1} + 1/n
    """
    # 递归实现(不很高效，但直接看逻辑最简单)
    if n == 0:
        return math.log((a + 1) / a)
    else:
        return -a * method1(a, n - 1) + 1 / n


def method2(a, n):
    """
    分段公式:
      若 a >= n/(n+1):
         I_n = (2a + 1)/( 2a(a+1)(n+1) )
      否则:
         I_n = 0.5( 1/((a+1)(n+1)) + 1/n )
    """
    if a >= n / (n + 1):
        return (2 * a + 1) / (2 * a * (a + 1) * (n + 1))
    else:
        return 0.5 * (1 / ((a + 1) * (n + 1)) + 1 / n)


# ============ 3. 对比不同 n 的误差并输出表格/作图 ============

def compare_methods(a, max_n=10):
    """
    对比 n=0..max_n 的 method1, method2 与 数值积分 的结果，
    返回包含误差的 DataFrame，并画出误差曲线。
    """
    results = {
        "n": [],
        "numeric_I": [],
        "method1": [],
        "method2": [],
        "err_m1": [],
        "err_m2": []
    }
    for n in range(max_n + 1):
        i_true = numeric_I(a, n)  # 数值积分参考值
        i_m1 = method1(a, n)
        i_m2 = method2(a, n)
        err1 = abs(i_m1 - i_true)
        err2 = abs(i_m2 - i_true)
        results["n"].append(n)
        results["numeric_I"].append(i_true)
        results["method1"].append(i_m1)
        results["method2"].append(i_m2)
        results["err_m1"].append(err1)
        results["err_m2"].append(err2)

    df = pd.DataFrame(results)
    # 打印或保存到 CSV
    # df.to_csv(f"compare_a_{a}.csv", index=False, encoding='utf-8')

    # 画出 method1, method2 的误差曲线
    plt.figure()
    plt.title(f"Absolute Error vs n, a={a}")
    plt.plot(df["n"], df["err_m1"], marker='o', label="Error of method1")
    plt.plot(df["n"], df["err_m2"], marker='s', label="Error of method2")
    plt.xlabel("n")
    plt.ylabel("Absolute Error")
    plt.yscale("log")  # 用对数坐标更容易看出差距
    plt.legend()
    plt.show()

    return df


# ============ 4. 分别测试 a=0.05 (小) 和 a=15 (大) ============

df_min = compare_methods(a=0.05, max_n=10)
print("===== a = 0.05 的对比结果 =====")
print(df_min)

df_max = compare_methods(a=15, max_n=10)
print("\n===== a = 15 的对比结果 =====")
print(df_max)
