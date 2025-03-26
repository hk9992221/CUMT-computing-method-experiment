import math

# 定义目标函数 f(x) = x^3 - sin(x) - 12*x + 1
def f(x):
    return x**3 - math.sin(x) - 12*x + 1

# f 的导数，用于牛顿法
def df(x):
    return 3*x**2 - math.cos(x) - 12

# 牛顿法实现（达到最大迭代次数时直接返回当前值），打印达到容忍度时的迭代次数
def newton_method(x0, tol=1e-6, max_iter=100):
    x = x0
    for i in range(1, max_iter+1):
        fx = f(x)
        dfx = df(x)
        if abs(dfx) < 1e-12:  # 避免除零
            print(f"牛顿法在第 {i} 步遇到导数过小，提前退出")
            return x
        x_new = x - fx / dfx
        if abs(x_new - x) < tol:
            print(f"牛顿法在第 {i} 步达到容忍度条件")
            return x_new
        x = x_new
    return x

# 计算立方根（注意负数）
def cube_root(val):
    if val >= 0:
        return val**(1/3)
    else:
        return -((-val)**(1/3))

# 不动点迭代法 - 方法 a: 使用 phi_1(x) = (12*x + sin(x) - 1)^(1/3)
def fixed_point_phi1(x0, tol=1e-6, max_iter=100, verbose=True):
    x = x0
    for i in range(1, max_iter+1):
        val = 12*x + math.sin(x) - 1
        x_new = cube_root(val)
        if abs(x_new - x) < tol:
            if verbose:
                print(f"phi_1 不动点法在第 {i} 步达到容忍度条件")
            return x_new
        x = x_new
    return x

# 不动点迭代法 - 方法 b: 使用 phi_2(x) = (x^3 - sin(x) + 1) / 12
def fixed_point_phi2(x0, tol=1e-6, max_iter=100, verbose=True):
    x = x0
    for i in range(1, max_iter+1):
        x_new = (x**3 - math.sin(x) + 1) / 12
        if abs(x_new - x) < tol:
            if verbose:
                print(f"phi_2 不动点法在第 {i} 步达到容忍度条件")
            return x_new
        x = x_new
    return x

# 埃特金加速：对不动点迭代函数进行加速迭代，内部调用 phi 时关闭提示信息
def aitken_acceleration(phi, x0, tol=1e-6, max_iter=100):
    x0_ait = x0
    for i in range(1, max_iter+1):
        # 关闭内部提示信息，避免重复打印
        x1 = phi(x0_ait, tol=tol, max_iter=max_iter, verbose=False)
        x2 = phi(x1, tol=tol, max_iter=max_iter, verbose=False)
        denominator = x2 - 2*x1 + x0_ait
        if abs(denominator) < 1e-12:
            print(f"埃特金加速法在第 {i} 步分母过小，提前退出")
            return x2
        x_acc = x0_ait - (x1 - x0_ait)**2 / denominator
        if abs(x_acc - x0_ait) < tol:
            print(f"埃特金加速法在第 {i} 步达到容忍度条件")
            return x_acc
        x0_ait = x_acc
    return x0_ait

# 主程序
if __name__ == "__main__":
    tol = 1e-6

    print("使用牛顿法求解：")
    # 取三个初值分别求根：[-3.5, 0.1, 3.5]
    x0_newton = [-3.5, 0.1, 3.5]
    for x0 in x0_newton:
        root = newton_method(x0, tol=tol)
        print(f"{root:12.8f}")

    print("\n使用 phi_1 不动点迭代法求解（适用于区间 [-4,-3] 与 [3,4]）：")
    x0_phi1 = [-3.5, 3.5]
    for x0 in x0_phi1:
        root = fixed_point_phi1(x0, tol=tol)
        print(f"{root:12.8f}")

    print("\n使用 phi_2 不动点迭代法求解（适用于区间 [0,0.2]）：")
    root_phi2 = fixed_point_phi2(0.1, tol=tol)
    print(f"{root_phi2:12.8f}")

    print("\n使用埃特金加速法对 phi_2 迭代序列求解（适用于区间 [0,0.2]）：")
    root_aitken = aitken_acceleration(fixed_point_phi2, 0.1, tol=tol)
    print(f"{root_aitken:12.8f}")
