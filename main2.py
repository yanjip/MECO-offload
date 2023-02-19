'''
@Author  ：Yan JP
@Created on Date：2023/2/12 16:54 
'''
import numpy as np
from UE import *
import define
from plot import *

def Picture1():
    res2=[]
    res3=[]
    aver=[]
    define.F_MEC=61e9
    define.UE_n=25

    for t in T:
        define.T_slot=t
        # u1 = UE_All()
        # u1.generate_ue()

        # -----------算法一
        u1.testAlgor1()
        s = u1.energy_all()
        # 算法2
        u1.lk_star2,u1.tk_star2=Alorithm2(u1.lk_star, u1.Ck, define.F_MEC,u1)
        print("-----------------------------------\n",u1.lk_star2,u1.tk_star2)
        print("---------------\n 算法二执行完毕")
        s2=u1.energy_all_2()

        u1.lk_star3, u1.tk_star3 = Alorithm3(u1.lk_star, u1.Ck, define.F_MEC, u1, define.T_slot)
        print("---------------\n 算法三执行完毕")

        s3 = u1.energy_all_3()
        print("算法三耗能：", s3)


        # cvxpy
        tk = [define.T_slot / u1.N] * u1.N
        s_equal = cvxSolve(u1.hk, u1.mk, u1.Rk, u1.Pk, u1.Ck, tk,define.F_MEC)

        res2.append(s2)
        res3.append(s3)
        aver.append(s_equal)
        print("算法二耗能：", s)
        print("equel耗能：", s_equal)

        print("t1,t2:", u1.get_F_bound())
        print("T_max:", u1.get_T_max())
        print("T_min:", get_Tmin(u1))

        with open('res.txt', 'a+') as F:
            F.write("T_slot:"+str(t)+"：\n算法二耗能：" + str(s2)+"      ")
            F.write("：算法三耗能：" + str(s3)+"      ")
            F.write("：equel耗能："+str(s_equal)+"      \n")
    plot_1(res2,res3,aver,T)

def Picture2():
    res2=[]
    res3=[]
    aver=[]
    define.F_MEC=82e9  #66还行
    define.T_slot=1.30

    for n in N:
        define.UE_n=n
        # u1 = UE_All()
        # u1.generate_ue()
        u1.N=n
        # -----------算法一
        u1.testAlgor1()
        s = u1.energy_all()
        # 算法2
        u1.lk_star2,u1.tk_star2=Alorithm2(u1.lk_star, u1.Ck, define.F_MEC,u1)
        print("-----------------------------------\n",u1.lk_star2,u1.tk_star2)
        print("---------------\n 算法二执行完毕")
        s2=u1.energy_all_2()

        u1.lk_star3, u1.tk_star3 = Alorithm3(u1.lk_star, u1.Ck, define.F_MEC, u1, define.T_slot)
        print("---------------\n 算法三执行完毕")

        s3 = u1.energy_all_3()
        print("算法三耗能：", s3)

        # cvxpy
        tk = [define.T_slot / n] * n
        s_equal = cvxSolve(u1.hk[:n], u1.mk[:n], u1.Rk[:n], u1.Pk[:n], u1.Ck[:n], tk[:n],define.F_MEC)

        res2.append(s2)
        res3.append(s3)
        aver.append(s_equal)


        with open('res.txt', 'a+') as F:
            F.write("UE_N:"+str(n)+"：\n算法二耗能：" + str(s2)+"      ")
            F.write("：算法三耗能：" + str(s3)+"      ")
            F.write("：equel耗能："+str(s_equal)+"      \n")

        print("t1,t2:", u1.get_F_bound())
        print("T_max:", u1.get_T_max())
        print("T_min:", get_Tmin(u1))
    plot_2(res2,res3,aver,N)
with open('res.txt', 'a+') as F:
    F.write("\n\n")

# T=np.array([0.55,0.65,0.75,0.85,0.95,1.5])+0.5
T=np.array([1.15,1.25,1.35,1.45,1.55])+0.0
# np.random.seed(6)
# u1 = UE_All()
# u1.generate_ue()
# Picture1()


# N=[10,15,20,25,30,35]
N=np.linspace(10,30,6,dtype=int)
np.random.seed(5) #5可以
define.UE_n=N[-1]
u1 = UE_All()
u1.generate_ue()
Picture2()

