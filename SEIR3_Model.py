import scipy.integrate as spi
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei']
# N为人群总数
N = 10000
# β为传染率系数
beta = 0.6
# gamma为恢复率系数
gamma = 0.1
# Te为疾病潜伏期
Te = 14
# I_0为感染者的初始人数
I_0 = 1
# E_0为潜伏者的初始人数
E_0 = 0
# R_0为治愈者的初始人数
R_0 = 0
# S_0为易感者的初始人数
S_0 = N - I_0 - E_0 - R_0
# T为传播时间
T = 150

# INI为初始状态下的数组
INI = (S_0, E_0, I_0, R_0)


def funcSEIR(inivalue, _):
    Y = np.zeros(4)
    X = inivalue
    # 易感个体变化
    Y[0] = - (beta * X[0] * X[2]) / N
    # 潜伏个体变化（每日有一部分转为感染者）
    Y[1] = (beta * X[0] * X[2]) / N - X[1] / Te
    # 感染个体变化
    Y[2] = X[1] / Te - gamma * X[2]
    # 治愈个体变化
    Y[3] = gamma * X[2]
    return Y


T_range = np.arange(0, T + 1)

RES = spi.odeint(funcSEIR, INI, T_range)

plt.plot(RES[:, 0], color='g', label='易感染者——Susceptible')
plt.plot(RES[:, 1], color='orange', label='潜伏者——Exposed')
plt.plot(RES[:, 2], color='r', label='传染者——Infection')
plt.plot(RES[:, 3], color='b', label='康复者——Recover')

plt.title('SEIR Model(具有潜伏期)')
plt.legend()
plt.xlabel('天数')
plt.ylabel('人数')
plt.show()
