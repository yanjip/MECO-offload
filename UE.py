'''
@Author  ：Yan JP
@Created on Date：2023/2/4 13:17 
'''
import define
import math
import numpy as np
from Algorithm import *
from cvx_solve import cvxSolve,solve_equel
class UE:
    def __init__(self,Pk,Rk,Fk,Ck,hk,id):
        self.id=id
        self.Pk=Pk
        self.Rk=Rk
        self.Fk=Fk
        self.Ck=Ck
        # self.lk=lk
        self.hk=hk
        self.mk=max(self.Rk-int(self.Fk*define.T_slot/self.Ck),0)
        self.vk=define.B*Ck*Pk*hk**2/(define.N0*np.log(2))
        print("vk:",self.vk)

        pass
    def request(self):
        pass
    def get_mk(self):
        return self.mk
    def phy(self):  #Priority Function
        if self.vk<1:
            self.phyk=0     #Priority Function
        else:
            self.phyk=define.N0*(self.vk*np.log(self.vk)-self.vk+1)/self.hk**2
        return self.phyk

    def get_new_phyk(self,u):
        self.vk2=define.B*self.Ck*(self.Pk-u)*self.hk**2/(define.N0*np.log(2))
        if self.vk2<1:
            self.phyk2=0     #Priority Function
        else:
            self.phyk2=define.N0*(self.vk2*np.log(self.vk2)-self.vk2+1)/self.hk**2
        return self.phyk2

    def offload(self):

        pass
    def local_consume(self,lk):
        return int(self.Rk-lk)*int(self.Ck)*(self.Pk)
        pass

    def offload_consume(self,lk,tk):
        if tk==0: return 0
        return tk*define.fx(lk/tk)/self.hk**2


    def compute_rate(self):
        pass

    def power_consum(self):
        pass


class UE_All:
    def __init__(self):
        # self.id=id
        self.N=define.UE_n
        # np.random.seed(0)
        self.Pk=np.random.uniform(0,20e-11, self.N)
        self.Rk=np.random.randint(define.Rk[0],define.Rk[1],self.N)*define.KB
        self.Fk=np.random.choice(define.Fk,self.N)
        self.Ck=np.random.randint(500,1500,self.N)
        # print(self.Pk,self.Rk,self.Fk,self.Ck)
        self.phyk=[0]*self.N
        # print(self.Fk[1])
        # print(type(self.Pk)

        #生成信道数据
        # H = np.random.rayleigh(scale=1, size= self.N)*1e-3
        # H = np.random.rayleigh(scale=1, size= self.N)
        H=np.sqrt(np.random.exponential(1e-6,self.N))
        # print(H)
        # while np.min(H)==0:
        #     # np.random.seed(0)
        #     H = np.sqrt(np.random.exponential(1e-6, self.N))

        # u, s, v = np.linalg.svd(H[0])
        self.hk= H
        # print("hk:",self.hk)
        self.mk=[0]*self.N

        #计算vk
    def generate_ue(self):
        self.ues=[UE]*self.N
        for i in range(self.N):
            self.ues[i]=UE(self.Pk[i],self.Rk[i],self.Fk[i],self.Ck[i],self.hk[i],i)
            self.phyk[i]=self.ues[i].phy()
            self.mk[i]=self.ues[i].get_mk()
        print("phyk:",self.phyk)

    def get_new_phy_all(self,u):
        self.phyk2=[0]*self.N
        for i in range(self.N):
            self.phyk2[i]=self.ues[i].get_new_phyk(u)
        return self.phyk2


    def testAlgor1(self):
        lamda_max=max(self.phyk)
        print("lamdamax:",lamda_max)
        self.lk_star,self.tk_star=Alorithm1(lamda_max,define.T_slot,self.hk,define.N0,define.B,self.phyk,self.mk,self.Rk)
        # print("------u1.lk_star, u1.tk_star:-----------------------------\n", u1.lk_star, u1.tk_star)

    # def testAlgor2(self):
    #     return Alorithm2(self.lk_star,self.Ck,define.F_MEC)

    def get_vk(self):
        self.vk=[0]*self.N
        for i in range(self.N):
            self.vk[i]=self.ues[i].vk
        pass
        return self.vk
    def energy_all(self):
        #本地耗能
        energy_local=[0]*self.N
        energy_offload=[0]*self.N

        for i in range(self.N):
            energy_local[i]=self.ues[i].local_consume(self.lk_star[i])
            energy_offload[i]=self.ues[i].offload_consume(self.lk_star[i],self.tk_star[i])
        self.energy_sum=sum(energy_local)+sum(energy_offload)
        return self.energy_sum

        pass
    # def sort(self):
    # def equal_allocation(self):
    def energy_all_equal_allocation(self,lk,tk):    #不用的
        # tk=[define.T_slot/self.N]*self.N
        #本地耗能
        energy_local=[0]*self.N
        energy_offload=[0]*self.N

        for i in range(self.N):
            energy_local[i]=self.ues[i].local_consume(self.lk[i])
            energy_offload[i]=self.ues[i].offload_consume(self.lk_star[i],self.t[i])
        self.energy_sum=sum(energy_local)+sum(energy_offload)
        return self.energy_sum

    def energy_all_2(self):
        #本地耗能
        energy_local=[0]*self.N
        energy_offload=[0]*self.N

        for i in range(self.N):
            energy_local[i]=self.ues[i].local_consume(self.lk_star2[i])
            energy_offload[i]=self.ues[i].offload_consume(self.lk_star2[i],self.tk_star2[i])
        self.energy_sum=sum(energy_local)+sum(energy_offload)
        return self.energy_sum

    def energy_all_3(self):
        #本地耗能
        energy_local=[0]*self.N
        energy_offload=[0]*self.N

        for i in range(self.N):
            energy_local[i]=self.ues[i].local_consume(self.lk_star3[i])
            energy_offload[i]=self.ues[i].offload_consume(self.lk_star3[i],self.tk_star3[i])
        self.energy_sum=sum(energy_local)+sum(energy_offload)
        return self.energy_sum
    def get_F_bound(self):
        t1= 0
        t2=0
        for i in range(self.N):
            t1+=int(self.Ck[i])*int(self.mk[i])
            t2+=int(self.Ck[i])*int(self.Rk[i])
        return t1,t2

    def get_T_max(self):
        res=0.0
        lamda=min(self.phyk)
        # if lamda==0:
        #     print("lamda==0")
        #     exit(0)

        for i in range(self.N):
            if self.phyk[i] >0:
                lk= self.Rk[i]
            elif self.phyk[i] ==0:
                lk = self.mk[i]
            # tk[i]  = np.log(2) * lk[i]  /((lambertw((lamda * hk[i]  ** 2 - N0) / (N0 * np.e)).real+1)*B)
            # tk[i]  = np.log(2) * lk[i]  /(B* (lambertw((lamda * hk[i]  ** 2 - N0) / (N0 * np.e)).real)+1) #旧版本
            res += np.log(2) * lk / (define.B * (lambertw((lamda * self.hk[i] ** 2 - define.N0) / (define.N0 * np.e)).real + 1))
        return res
        pass
if __name__ == '__main__':
    np.random.seed(700)
    u1=UE_All()
    u1.generate_ue()

    #-----------算法一
    u1.testAlgor1()
    s=u1.energy_all()



    #算法2
    u1.lk_star2,u1.tk_star2=Alorithm2(u1.lk_star, u1.Ck, define.F_MEC,u1)
    print("-----------------------------------\n",u1.lk_star2,u1.tk_star2)
    print("---------------\n 算法二执行完毕")
    s2=u1.energy_all_2()
    print("算法二耗能：",s2)

    # # print("算法一耗能：",s)
    # print("equel耗能：",s_equal)
    # print("T_bound",get_Tmin(u1))
    #
    # #算法三
    u1.lk_star3,u1.tk_star3=Alorithm3(u1.lk_star, u1.Ck, define.F_MEC,u1,define.T_slot)
    print("---------------\n 算法三执行完毕")

    # # print("-----------------------------------\n",tk,lk)
    s3=u1.energy_all_3()
    print("算法三耗能：",s3)
    print("算法二耗能：",s2)
    print("算法一耗能：", s)

    #cvxpy

    tk=[define.T_slot/u1.N]*u1.N
    s_equal=cvxSolve(u1.hk,u1.mk,u1.Rk,u1.Pk,u1.Ck,tk,define.F_MEC) #函数里面打印了结果


    #equel_new
    # s_equal2=solve_equel(u1)
    # print("equel耗能：",s_equal2)


    print("t1,t2:",u1.get_F_bound())
    print("T_max:",u1.get_T_max())
    print("T_min:",get_Tmin(u1))


    #---------------------------------
    # u1.lk_star=u1.mk
    # _,u1.tk_star=compute_tk()
    # s=u1.energy_all()
    # print("耗能：",s)
