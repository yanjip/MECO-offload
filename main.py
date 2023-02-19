'''
@Author  ：Yan JP
@Created on Date：2023/2/4 13:19 
'''
import numpy as np
from UE import *
import define
from plot import *


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
def Picture3():
    alor_3=[]
    aver=[]

    for f_mec in F:
        define.F_MEC=f_mec
        # np.random.seed(3)
        # u1 = UE_All()
        # u1.generate_ue()
        # # -----------算法一
        # u1.testAlgor1()
        # s = u1.energy_all()

        u1.lk_star3, u1.tk_star3 = Alorithm3(u1.lk_star, u1.Ck, f_mec, u1, define.T_slot)
        print("---------------\n 算法三执行完毕")

        # # print("-----------------------------------\n",tk,lk)
        s3 = u1.energy_all_3()
        print("算法三耗能：", s3)
        # print("算法二耗能：",s2)
        print("算法一耗能：", s)

        # cvxpy

        tk = [define.T_slot / u1.N] * u1.N
        s_equal = cvxSolve(u1.hk, u1.mk, u1.Rk, u1.Pk, u1.Ck, tk, f_mec)  # 函数里面打印了结果

        alor_3.append(s3)
        aver.append(s_equal)
    print("t1,t2:",u1.get_F_bound())
    print("T_max:",u1.get_T_max())
    print("T_min:",get_Tmin(u1))
    print(alor_3)
    print(aver)
    plot_3(alor_3,aver,F)


# with open('res.txt', 'a+') as F:
#     F.write("\n\n")
# Picture1()
#
# with open('res.txt', 'a+') as F:
#     F.write("\n\n")
# Picture2()
# T=[0.5,0.6,0.7,0.8,0.9,1.0]
# N=[10,15,20,25,30,35]
F=[56e9,58e9,60e9,62e9,64e9,66e9,68e9,72e9,76e9,80e9]
F=np.array(F)-10e9  #原本18
define.T_slot=1.05
define.UE_n=25
np.random.seed(700)
u1 = UE_All()
u1.generate_ue()
# -----------算法一
u1.testAlgor1()
s = u1.energy_all()

Picture3()

# u1.lk_star2,u1.tk_star2=Alorithm2(u1.lk_star, u1.Ck, define.F_MEC,u1)
# print("-----------------------------------\n",u1.lk_star2,u1.tk_star2)

# print("---------------\n 算法二执行完毕")
# s2=u1.energy_all_2()
# print("算法二耗能：",s2)
# print("t1,t2:", u1.get_F_bound())
# print("T_max:", u1.get_T_max())
# print("T_min:", get_Tmin(u1))

