import scipy.integrate as spi
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']
# N为人群总数
N = 10000
# β为传染率系数
B = 0.05
# y为恢复健康的概率
y = 0.1

# R为康复者人数
R = 0
# I为感染者的初始人数
I = 1
# S为易感者的初始人数
S = N - I - R
# 感染者每天接触人数
r = 10
# T为传播时间
T = 100

# INI为初始状态下的数组
INI = (S, I, R)


def funcSIR(inivalue, _):
    Y = np.zeros(3)
    X = inivalue
    # 易感个体S变化
    Y[0] = - (r * B * X[0] * X[1]) / N
    # 感染个体I变化
    Y[1] = (r * B * X[0] * X[1]) / N - y * X[1]
    # 康复个体R变化
    Y[2] = y * X[1]
    return Y


T_range = np.arange(0, T + 1)

RES = spi.odeint(funcSIR, INI, T_range)
# odeint()函数需要至少三个变量，第一个是微分方程函数，第二个是微分方程初值，第三个是微分的自变量。

plt.plot(RES[:, 0], color='g', label='易感染者——Susceptible')
plt.plot(RES[:, 1], color='r', label='传染者——Infection')
plt.plot(RES[:, 2], color='b', label='康复者——Recover')
plt.title('SIR Model')
plt.legend()
plt.xlabel('天数')
plt.ylabel('人数')
plt.show()
