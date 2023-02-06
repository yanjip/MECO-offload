'''
@Author  ：Yan JP
@Created on Date：2023/2/5 23:43 
'''
import cvxpy as cvx
from define import *
def cvxSolve(hk,mk,Rk,Pk,Ck):
    lk=cvx.Variable(UE_n)
    tk=cvx.Variable(UE_n)

    question=cvx.sum(
        cvx.multiply(1/hk**2,cvx.multiply(tk,
                                          (cvx.exp(cvx.multiply(np.log(2),cvx.multiply(lk,1/tk)/B)-1))*N0
                                          )
    )+cvx.multiply(Pk,cvx.multiply(Ck,(Rk-lk)))
    )
    obj=cvx.Minimize(question)
    constrains=[cvx.sum(tk)<=T_slot,tk>=0.0,Rk>=lk,lk>=mk]
    prob=cvx.Problem(obj,constrains)
    ans=prob.solve(solver=cvx.MOSEK)

    print(prob.value)