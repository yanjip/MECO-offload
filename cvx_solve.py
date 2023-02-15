'''
@Author  ：Yan JP
@Created on Date：2023/2/5 23:43 
'''
import cvxpy as cvx
# from define import *
import define
import numpy as np
# from UE import  UE_All
def cvxSolve(hk,mk,Rk,Pk,Ck,tk,F):
    tk=np.array(tk)
    lk=cvx.Variable(define.UE_n)
    # tk=cvx.Variable(UE_n)


    question=cvx.sum(
        cvx.multiply(1/hk**2,
                     cvx.multiply(tk,
                                      (cvx.exp(cvx.multiply(np.log(2),
                                                            cvx.multiply(lk,1/tk)/define.B))-1)*define.N0
                                      )
    )+cvx.multiply(Pk,cvx.multiply(Ck,(Rk-lk)))
    )
    obj=cvx.Minimize(question)
    # constrains=[cvx.sum(tk)<=T_slot,tk>=0.0,Rk>=lk,lk>=mk]
    # constrains=[Rk>=lk,lk>=mk]

    constrains=[Rk>=lk,lk>=mk,cvx.sum(cvx.multiply(Ck,lk))<=F]

    prob=cvx.Problem(obj,constrains)
    # prob.solve(solver=cvx.MOSEK)
    prob.solve()


    print("equel allocation:",prob.value)
    # print("lk:",lk.value)
    return prob.value

def solve_equel(u1):
    tk=[define.T_slot/u1.N]*u1.N
    lk=u1.Rk
    vk=u1.get_vk()
    F_temp=0.0
    for i in range(u1.N):
        F_temp+=int(u1.Ck[i])*int(u1.Rk[i])

        if vk[i]<1 or F_temp>define.F_MEC: #卸载mk
            # tk[i]=0
            lk[i]=u1.mk[i]
            F_temp -= int(u1.Ck[i]) * int(u1.Rk[i])
            F_temp += int(u1.Ck[i])*int(u1.mk[i])

        # else:  #可以卸载Rk
        #     lk[i]=u1.Rk[i]
    # print("F_TEMP",F_temp)
    # t1= 0
    # t2=0
    # for i in range(u1.N):
    #     t1+=int(u1.Ck[i])*int(u1.mk[i])
    #     t2+=int(u1.Ck[i])*int(u1.Rk[i])
            # F_temp+=int(u1.Ck[i])*int(lk[i])
            # if F_temp>define.F_MEC:  #不能卸载，lk取最大
            #     F_temp -= int(u1.Ck[i]) * int(lk[i])
                # tk[i]=0
            # else:
            #     lk[i]=u1.mk[i]

    energy=0.0
    # print(lk)
    # print(tk)
    for i in range(u1.N):
        E_loc=(u1.Rk[i]-lk[i])*u1.Ck[i]*u1.Pk[i]
        if tk[i]!=0:
            E_off=tk[i]*define.fx(lk[i]/tk[i])/u1.hk[i]**2
        else:
            E_off=0
        energy+=(E_off+E_loc)

    return energy


    pass