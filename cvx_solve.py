'''
@Author  ：Yan JP
@Created on Date：2023/2/5 23:43 
'''
import cvxpy as cvx
# from define import *
import define
import numpy as np
def cvxSolve(hk,mk,Rk,Pk,Ck,tk):
    tk=np.array(tk)
    lk=cvx.Variable(define.UE_n)
    # tk=cvx.Variable(UE_n)


    question=cvx.sum(
        cvx.multiply(1/hk**2,cvx.multiply(tk,
                                          (cvx.exp(cvx.multiply(np.log(2),cvx.multiply(lk,1/tk)/define.B))-1)*define.N0
                                          )
    )+cvx.multiply(Pk,cvx.multiply(Ck,(Rk-lk)))
    )
    obj=cvx.Minimize(question)
    # constrains=[cvx.sum(tk)<=T_slot,tk>=0.0,Rk>=lk,lk>=mk]
    constrains=[Rk>=lk,lk>=mk]

    # constrains=[Rk>=lk,lk>=mk,cvx.sum(cvx.multiply(Ck,lk))<=F_MEC]

    prob=cvx.Problem(obj,constrains)
    prob.solve(solver=cvx.MOSEK)

    print("equel allocation:",prob.value)
    # print("lk:",lk.value)
    return prob.value