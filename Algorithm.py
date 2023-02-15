'''
@Author  ：Yan JP
@Created on Date：2023/2/4 13:12 
'''
import numpy as np
from scipy.special import lambertw
import cvxpy as cvx
import time
from cvx_solve import cvxSolve
from UE import *

def compute_lk_tk(lamda,hk,N0,B,phyk,mk,Rk):
    lk=[0]*len(Rk)
    tk=[0]*len(Rk)

    for i in range(len(Rk)):
        if mk[i]==0:
            continue
        if phyk[i] < lamda:
            lk[i]  = mk[i]
        elif phyk[i]  > lamda:
            lk[i]  = Rk[i]
        # tk[i]  = np.log(2) * lk[i]  /((lambertw((lamda * hk[i]  ** 2 - N0) / (N0 * np.e)).real+1)*B)
        # tk[i]  = np.log(2) * lk[i]  /(B* (lambertw((lamda * hk[i]  ** 2 - N0) / (N0 * np.e)).real)+1) #旧版本
        tk[i] = np.log(2) * lk[i] / (B * (lambertw((lamda * hk[i] ** 2 - N0) / (N0 * np.e)).real+1))
    # print("*********************\n",tk)
    return sum(tk),lk,tk


def get_Tmin(u1):
    phy_max = max(u1.phyk)
    max_k = u1.phyk.index(phy_max)
    R=u1.Rk[max_k]
    C=u1.Ck[max_k]
    P=u1.Pk[max_k]
    h=u1.hk[max_k]
    return  R/(define.B*np.log2(define.B*C*P*h**2/(define.N0*np.log(2))))

def Alorithm1(lamda,T,hk,N0,B,phyk,mk,Rk):

    print("lamda_max_alor1:",lamda)


    lamda_l=0.0
    lamda_h=lamda
    T_l,_,_=compute_lk_tk(lamda_l,hk,N0,B,phyk,mk,Rk)
    T_h,_,_=compute_lk_tk(lamda_h,hk,N0,B,phyk,mk,Rk)
    print("T_l:", T_l) #nan  lambertw(-0.36787944117144235) 应该是小数点太多了，计算出来为nan
    print("T_h:", T_h)
    delta = 0.000000001
    if abs(T_l - T) < delta:
        lamda = lamda_l
        _, lk_star, tk_star = compute_lk_tk(lamda, hk, N0, B, phyk, mk, Rk)
        return lk_star, tk_star
    if abs(T_h - T) < delta:
        lamda = lamda_h
        _,lk_star,tk_star=compute_lk_tk(lamda,hk,N0,B,phyk,mk,Rk)
        return lk_star,tk_star

    while T_l!=T and T_h!=T:
        lamda_m = (lamda_l + lamda_h) / 2
        print("lamda_m:",lamda_m)
        # if lamda_m==lamda_l:
        #     return
        T_m,_,_ = compute_lk_tk(lamda_m, hk, N0, B, phyk, mk, Rk)
        print('T_m:',T_m)
        if abs(T_m-T)<delta:
            lamda=lamda_m
            break
        if T_m<T:
            lamda_h=lamda_m
        else:
            lamda_l=lamda_m
            print("lamda_l:", lamda_l)
            print("lamda_h:",lamda_h)

    _,lk_star,tk_star=compute_lk_tk(lamda,hk,N0,B,phyk,mk,Rk)
    return lk_star,tk_star

def compute_F(Ck,lk):
    res=0.0
    lk=np.array(lk)     #不写这个会报错：RuntimeWarning: overflow encountered in long_scalars
    for i in range(len(Ck)):
        res+=int(Ck[i])*int(lk[i])
    return res

    pass
def get_u_max(Pk,N0,B,Ck,hk):
    u_list=[0]*len(Ck)
    for i in range(len(Ck)):
        u_list[i]=Pk[i]-(N0*np.log(2))/(B*Ck[i]*hk[i]**2)
    return max(u_list)

def Alorithm2(lk_base1,Ck,F_MEC,ueall):
    F=compute_F(Ck,lk_base1)
    print("F:",F)
    if F<=F_MEC:
        return lk_base1,ueall.tk_star
    u_l=0
    u_h=get_u_max(ueall.Pk,define.N0,define.B,ueall.Ck,ueall.hk)

    ueall.get_new_phy_all(u_l)  #每个ue计算了自己的pyhk2
    lamda_max = max(ueall.phyk2)
    lk_algor1_ul, tk_star = Alorithm1(lamda_max, define.T_slot, ueall.hk, define.N0, define.B, ueall.phyk2, ueall.mk,ueall.Rk)
    F_l=compute_F(Ck,lk_algor1_ul)


    # ueall.get_new_phy_all(u_h)  #每个ue计算了自己的pyhk2
    # lamda_max = max(ueall.phyk2)
    # print("lamdamax:", lamda_max)
    # lk_algor1_uh, tk_star = Alorithm1(lamda_max, define.T_slot, ueall.hk, define.N0, define.B, ueall.phyk2, ueall.mk,ueall.Rk)
    # F_h=compute_F(Ck,lk_algor1_uh)

    while F_l!=F_MEC :
        u_m=(u_l+u_h)/2
        ueall.get_new_phy_all(u_m)  # 每个ue计算了自己的pyhk2
        lamda_max = max(ueall.phyk2)
        lk_algor1_um, tk_star = Alorithm1(lamda_max, define.T_slot, ueall.hk, define.N0, define.B, ueall.phyk2,
                                          ueall.mk, ueall.Rk)
        F_m=compute_F(Ck,lk_algor1_um)
        delta=0.1e9
        print("F_M:",F_m)
        print("差值：",abs(F_MEC-F_m))
        if abs(F_MEC-F_m)<delta:
            return lk_algor1_um, tk_star
        if F_m<F_MEC:
            u_h=u_m
        else:
            u_l=u_m



    # return lk_star, tk_star


    pass
#-------------------baseline---------------------#
# def Equal_allocation(u1:UE_All):

def compute_tk(lamda,lk,B,hk,N0):
    tk=[0]*len(lk)
    for i in range(len(tk)):
        tk[i] = np.log(2) * lk[i] / (B * (lambertw((lamda * hk[i] ** 2 - N0) / (N0 * np.e)).real+1))
    return sum(tk),tk

def Alorithm3(lk_base1,Ck,F_MEC,ueall,T):
    lk=[0]*len(Ck)
    tk=[0]*len(Ck)
    F=compute_F(Ck,lk_base1)
    print("F:",F)
    F_min=compute_F(Ck,ueall.mk)
    print("F_min:",F_min)
    if F<=F_MEC:
        return lk_base1,ueall.tk_star

    idx=[i for i in range(ueall.N)]
    d=dict(zip(idx,ueall.phyk))
    d=sorted(d.items(), key=lambda x: x[1], reverse=True) #现在是列表

    F_temp=0
    for i in d:
        F_temp += int(Ck[i[0]]) * int(ueall.Rk[i[0]])
        if F_temp<=F_MEC:
            lk[i[0]]=ueall.Rk[i[0]]
        else:
            F_temp -= int(Ck[i[0]]) * int(ueall.Rk[i[0]])
            lk[i[0]]=ueall.mk[i[0]]
            F_temp += int(Ck[i[0]]) * int(ueall.mk[i[0]])
    print("F_temp:",F_temp)
    lamda_max=max(ueall.phyk)

    # Alorithm1(lamda_max,define.T_slot,ueall.hk,define.N0,define.B,
    #           ueall.phyk,ueall.mk,ueall.Rk)

    lamda_l = 0.0
    lamda_h = lamda_max
    T_l, _ = compute_tk(lamda_l,lk,define.B,ueall.hk,define.N0)
    T_h, _= compute_tk(lamda_h,lk,define.B,ueall.hk,define.N0)
    print("T_l:", T_l)  # nan  lambertw(-0.36787944117144235) 应该是小数点太多了，计算出来为nan
    print("T_h:", T_h)
    delta = 0.0000001
    while T_l != T and T_h != T:
        lamda_m = (lamda_l + lamda_h) / 2
        print("lamda_m:", lamda_m)
        # if lamda_m==lamda_l:
        #     return
        T_m, tk = compute_tk(lamda_m, lk, define.B, ueall.hk, define.N0)

        print('T_m:', T_m)
        if abs(T_m - T) < delta:
            lamda = lamda_m
            break
        if T_m < T:
            lamda_h = lamda_m
        else:
            lamda_l = lamda_m
            print("lamda_l:", lamda_l)
            print("lamda_h:", lamda_h)
    # T_m, _ = compute_tk(lamda_m, lk, define.B, ueall.hk, define.N0)

    # _, lk_star, tk_star = compute_lk_tk(lamda, hk, N0, B, phyk, mk, Rk)
    return  lk,tk



# if __name__ == '__main__':
