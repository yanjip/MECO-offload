'''
@Author  ：Yan JP
@Created on Date：2023/2/4 13:17 
'''
import define
import math
import numpy as np
from Algorithm import *
from cvx_solve import cvxSolve
class UE:
    def __init__(self,Pk,Rk,Fk,Ck,hk,id):
        self.id=id
        self.Pk=Pk
        self.Rk=Rk
        self.Fk=Fk
        self.Ck=Ck
        # self.lk=lk
        self.hk=hk
        self.mk=max(self.Rk-self.Fk*define.T_slot/self.Ck,0)
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
            self.phyk=define.N0*(self.vk*np.log(self.vk)-self.vk+1)
        return self.phyk

    def offload(self):

        pass
    def local_consume(self,lk):
        return (self.Rk-lk)*self.Ck*self.Pk
        pass

    def offload_consume(self,lk,tk):
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
        self.Rk=np.random.randint(define.Rk[0],define.Rk[1],self.N)
        self.Fk=np.random.choice(define.Fk,self.N)
        self.Ck=np.random.randint(500,1500,self.N)
        # print(self.Pk,self.Rk,self.Fk,self.Ck)
        self.phyk=[0]*self.N
        # print(self.Fk[1])
        # print(type(self.Pk)

        #生成信道数据
        H = np.random.rayleigh(scale=2, size= self.N)*1e-3
        # H = np.random.rayleigh(scale=2, size= self.N)
        # print('标准差为2瑞利分布：\n', H)
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

    def testAlgor1(self):
        lamda_max=max(self.phyk)
        print("lamdamax:",lamda_max)
        self.lk_star,self.tk_star=Alorithm1(lamda_max,define.T_slot,self.hk,define.N0,define.B,self.phyk,self.mk,self.Rk)


    def offload_user(self):

        pass
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
    def energy_all_equal_allocation(self,lk,tk):
        # tk=[define.T_slot/self.N]*self.N
        #本地耗能
        energy_local=[0]*self.N
        energy_offload=[0]*self.N

        for i in range(self.N):
            energy_local[i]=self.ues[i].local_consume(self.lk[i])
            energy_offload[i]=self.ues[i].offload_consume(self.lk_star[i],self.t[i])
        self.energy_sum=sum(energy_local)+sum(energy_offload)
        return self.energy_sum

if __name__ == '__main__':
    np.random.seed(0)
    u1=UE_All()
    u1.generate_ue()
    u1.testAlgor1()
    s=u1.energy_all()
    print("耗能：")
    print(s)

    #cvxpy
    cvxSolve(u1.hk,u1.mk,u1.Rk,u1.Pk,u1.Ck)

