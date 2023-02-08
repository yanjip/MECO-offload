'''
@Author  ：Yan JP
@Created on Date：2023/2/4 13:19 
'''
import numpy as np
from UE import *
import define

T=[0.5,0.6,0.7,0.8,0.9,1.0]
N=[10,15,20,25,30,35]
def Picture1():
    np.random.seed(23)

    for t in T:
        define.T_slot=t
        u1 = UE_All()
        u1.generate_ue()

        # -----------算法一
        u1.testAlgor1()
        s = u1.energy_all()

        # cvxpy
        tk = [define.T_slot / u1.N] * u1.N
        s_equal = cvxSolve(u1.hk, u1.mk, u1.Rk, u1.Pk, u1.Ck, tk)

        print("算法一耗能：", s)
        print("equel耗能：", s_equal)

        with open('res.txt', 'a+') as F:
            F.write("T_slot:"+str(t)+"：算法一耗能：" + str(s)+"      ")
            F.write("T_slot:"+str(t)+"：equel耗能："+str(s_equal)+"      \n")

def Picture2():
    np.random.seed(23)

    for n in N:
        define.UE_n=n
        u1 = UE_All()
        u1.generate_ue()

        # -----------算法一
        u1.testAlgor1()
        s = u1.energy_all()

        # cvxpy
        tk = [define.T_slot / u1.N] * u1.N
        s_equal = cvxSolve(u1.hk, u1.mk, u1.Rk, u1.Pk, u1.Ck, tk)

        with open('res.txt', 'a+') as F:
            F.write("UE_N:"+str(n)+"：算法一耗能：" + str(s)+"      ")
            F.write("UE_N:"+str(n)+"：equel耗能："+str(s_equal)+"      \n")

with open('res.txt', 'a+') as F:
    F.write("\n\n")
Picture1()

with open('res.txt', 'a+') as F:
    F.write("\n\n")
Picture2()
# u1 = UE_All()
# u1.generate_ue()
#
# # -----------算法一
# u1.testAlgor1()
# s = u1.energy_all()
#
# # cvxpy
# tk = [define.T_slot / u1.N] * u1.N
# s_equal = cvxSolve(u1.hk, u1.mk, u1.Rk, u1.Pk, u1.Ck, tk)
#
# print("算法一耗能：", s)
# print("equel耗能：", s_equal)