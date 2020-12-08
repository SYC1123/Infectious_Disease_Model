import scipy.integrate as spi
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']
# β为传染率系数
B = 0.03
# y为恢复健康的概率
y = 0.1
# a为潜伏者转化率
a = 0.1
# r2潜伏者接触人数
r2 = 20
# B2为潜伏者传染率
B2 = 0.03

# N为人群总数
N = 10000
# E为潜伏者人数
E = 0
# R为康复者人数
R = 0
# I为感染者的初始人数
I = 1
# S为易感者的初始人数
S = N - I - E - R
# 感染者每天接触人数
r = 20
# T为传播时间
T = 150

# INI为初始状态下的数组
INI = (S, E, I, R)


def funcSEIR(inivalue, _):
    Y = np.zeros(4)
    X = inivalue
    # 易感个体S变化
    Y[0] = - (r * B * X[0] * X[2]) / N - (r2 * B2 * X[1] * X[0]) / N
    # 潜伏者个体E变化
    Y[1] = (r * B * X[0] * X[2]) / N - a * X[1] + (r2 * B2 * X[1] * X[0]) / N
    # 感染个体I变化
    Y[2] = a * X[1] - y * X[2]
    # 康复个体R变化
    Y[3] = y * X[2]

    return Y


T_range = np.arange(0, T + 1)

RES = spi.odeint(funcSEIR, INI, T_range)
# odeint()函数需要至少三个变量，第一个是微分方程函数，第二个是微分方程初值，第三个是微分的自变量。

plt.plot(RES[:, 0], color='g', label='易感染者——Susceptible')
plt.plot(RES[:, 1], color='orange', label='潜伏者——Exposed')
plt.plot(RES[:, 2], color='r', label='传染者——Infection')
plt.plot(RES[:, 3], color='b', label='康复者——Recover')
plt.title('SEIR Model(Exposed可传染)')
plt.legend()
plt.xlabel('天数')
plt.ylabel('人数')
plt.show()
