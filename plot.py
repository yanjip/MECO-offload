'''
@Author  ：Yan JP
@Created on Date：2023/2/4 13:06 
'''
import matplotlib.pyplot as plt
import scipy.io as sio                     # import scipy.io for .mat file I/
import numpy as np

def plot_3(alor3,aver,F):
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.plot(F,alor3,linestyle='--',label='算法三')
    plt.plot(F,aver,linestyle='--',label='baseline')
    plt.legend(loc='best')
    plt.xlabel('Cloud computaion capacity')
    plt.ylabel('total mobile energy consumption(J)')
    plt.show()

def plot_2(res2,res3,aver,N):
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.plot(N,res2,linestyle='--',label='算法三')
    plt.plot(N,res3,label='算法二')
    plt.plot(N,aver,linestyle='--',label='baseline')
    plt.legend(loc='best')
    plt.xlabel('Time Slot duration(S)')
    plt.ylabel('total mobile energy consumption(J)')
    plt.show()

def plot_1(res2,res3,aver,T):
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.plot(T,res2,linestyle='--',label='算法三')
    plt.plot(T,res3,linestyle='--',label='算法二')
    plt.plot(T,aver,linestyle='--',label='baseline')
    plt.legend(loc='best')
    plt.xlabel('Time Slot duration(S)')
    plt.ylabel('total mobile energy consumption(J)')
    plt.show()


if __name__ == '__main__':
    a=[0.2329720225338956, 0.22596623756503365, 0.1764518870890157, 0.17192083017452436, 0.1684701227822357,
     0.1684701227822357, 0.1684701227822357]
    b=[1.8562316046094558, 1.315696218494044, 0.9603252402803788, 0.8098304181433332, 0.8040087946419409,
     0.8040087424254662, 0.8040087208885199]
    F = [56e9, 60e9, 64e9, 68e9, 72e9, 76e9, 80e9]

    plot_3(a,b,F)

# plt.ylim(0,2000000)
# plt.plot(x, c, marker='s',label='最优解',markerfacecolor='none',clip_on=False,linewidth = 1, ms = 7)
# plt.plot(x, c2, marker='o',label='深度强化学习',markerfacecolor='none',clip_on=False,linewidth = 1, ms = 7)
# plt.plot(x, c1, marker='d',label='完全卸载',markerfacecolor='none',clip_on=False,linewidth = 1, ms = 7,color = 'r')
# plt.plot(x, d, marker='*',label='完全本地计算',markerfacecolor='none',clip_on=False,linewidth = 1, ms = 7)