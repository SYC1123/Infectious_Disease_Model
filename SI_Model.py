import scipy.integrate as spi
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei']
# N为人群总数
N = 10000
# β为传染率系数
B = 0.01
# I为感染者的初始人数
I = 1
# S为易感者的初始人数
S = N - I
# 感染者每天接触人数
r = 10
# T为传播时间
T = 200

# INI为初始状态下的数组
INI = (S, I)


def funcSI(inivalue, _):
    Y = np.zeros(2)
    X = inivalue
    # 易感个体变化
    Y[0] = - (r * B * X[0] * X[1]) / N
    # 感染个体变化
    Y[1] = (r * B * X[0] * X[1]) / N
    return Y


T_range = np.arange(0, T + 1)

RES = spi.odeint(funcSI, INI, T_range)
# odeint()函数需要至少三个变量，第一个是微分方程函数，第二个是微分方程初值，第三个是微分的自变量。

plt.plot(RES[:, 0], color='g', label='易感染者——Susceptible')
plt.plot(RES[:, 1], color='r', label='传染者——Infection')
plt.title('SI Model')
plt.legend()
plt.xlabel('天数')
plt.ylabel('人数')
plt.show()
