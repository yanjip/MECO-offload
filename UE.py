'''
@Author  ：Yan JP
@Created on Date：2023/2/4 13:17 
'''
import define
import math
import numpy as np
class UE:
    def __init__(self,Pk,Rk,Fk,Ck,lk,hk):

        self.Pk=Pk
        self.Rk=Rk
        self.Fk=Fk
        self.Ck=Ck
        self.lk=lk
        self.hk=hk
        self.mk=max(self.Rk-self.Fk*define.T_slot/self.Ck,0)
        self.vk=define.B*Ck*Pk*hk**2/(define.N0*np.log(2))
        if self.vk<1:
            self.phyk=0     #Priority Function
        else:
            self.phyk=define.N0*(self.vk*np.log(self.vk)-self.vk+1)

        pass
    def request(self):
        pass
    def phy(self):  #Priority Function
        return self.phyk
    def offload(self):

        pass
    def local_consume(self):
        pass
    def offload_consume(self,tk,hk,lk):
        return (tk/hk**2)*define.fx(x=lk/tk)

    def compute_rate(self):
        pass

    def power_consum(self):
        pass


class UE_All:
    def __init__(self):
        self.N=define.UE_n
        np.random.seed(0)
        self.Pk=np.random.uniform(0,20e-11, self.N)
        self.Rk=np.random.randint(define.Rk[0],define.Rk[1],self.N)
        self.Fk=np.random.choice(define.Fk,self.N)
        self.Ck=np.random.randint(500,1500,self.N)
        print(self.Pk,self.Rk,self.Fk,self.Ck)
        # print(type(self.Pk))
    def offload_user(self):

        pass
    def energy_all(self):
        pass
if __name__ == '__main__':
    u1=UE_All()
