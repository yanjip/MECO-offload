'''
@Author  ：Yan JP
@Created on Date：2023/2/4 13:12 
'''
import numpy as np
from scipy.special import lambertw
import cvxpy as cvx
import time

def compute_lk_tk(lamda,hk,N0,B,phyk,mk,Rk):
    lk=[0]*len(Rk)
    tk=[0]*len(Rk)

    for i in range(len(Rk)):
        if phyk[i] < lamda:
            lk[i]  = mk[i]
        elif phyk[i]  > lamda:
            lk[i]  = Rk[i]
        tk[i]  = np.log(2) * lk[i]  /(B* (lambertw((lamda * hk[i]  ** 2 - N0) / (N0 * np.e)).real)+1)


    return sum(tk),lk,tk


def Alorithm1(lamda,T,hk,N0,B,phyk,mk,Rk):
    lamda_l=0
    lamda_h=lamda
    T_l,_,_=compute_lk_tk(lamda_l,hk,N0,B,phyk,mk,Rk)
    T_h,_,_=compute_lk_tk(lamda_h,hk,N0,B,phyk,mk,Rk)
    print("T_l:", T_l)
    print("T_h:", T_h)
    delta = 0.000001
    while T_l!=T and T_h!=T:
        lamda_m=(lamda_l+lamda_h)/2
        print("lamda_m:",lamda_m)
        T_m,_,_ = compute_lk_tk(lamda_m, hk, N0, B, phyk, mk, Rk)
        print('T_m:',T_m)
        if abs(T_m-T)<delta:
            lamda=lamda_m
            break
        if T_m<T:
            lamda_h=lamda_m
        else:
            lamda_l=lamda_m
    _,lk_star,tk_star=compute_lk_tk(lamda,hk,N0,B,phyk,mk,Rk)
    return lk_star,tk_star

#-------------------baseline---------------------#
def Equal_allocation(T,N,Rk,mk):


    pass

def optimal():

    pass

def sub_optimal():

    pass






def bisection(h, M, weights=[]):
    # the bisection algorithm proposed by Suzhi BI
    # average time to find the optimal: 0.012535839796066284 s

    # parameters and equations
    o = 100
    p = 3
    u = 0.51
    eta1 = ((u * p) ** (1.0 / 3)) / o
    ki = 10 ** -26
    eta2 = u * p / 10 ** -10
    B = 2 * 10 ** 6
    Vu = 1.1
    epsilon = B / (Vu * np.log(2))
    x = []  # a =x[0], and tau_j = a[1:]

    M0 = np.where(M == 0)[0]
    M1 = np.where(M == 1)[0]

    hi = np.array([h[i] for i in M0])
    hj = np.array([h[i] for i in M1])

    if len(weights) == 0:
        # default weights [1, 1.5, 1, 1.5, 1, 1.5, ...]
        weights = [1 if i % 2 == 1 else 1 for i in range(len(M))]

    wi = np.array([weights[M0[i]] for i in range(len(M0))])
    wj = np.array([weights[M1[i]] for i in range(len(M1))])



    def phi(v, j):
        return 1 / (-1 - 1 / (lambertw(-1 / (np.exp(1 + v / wj[j] / epsilon))).real))

    def p1(v):
        p1 = 0
        for j in range(len(M1)):
            p1 += hj[j] ** 2 * phi(v, j)

        return 1 / (1 + p1 * eta2)

    def Q(v):
        sum1 = sum(wi * eta1 * (hi / ki) ** (1.0 / 3)) * p1(v) ** (-2 / 3) / 3
        sum2 = 0
        for j in range(len(M1)):
            sum2 += wj[j] * hj[j] ** 2 / (1 + 1 / phi(v, j))
        return sum1 + sum2 * epsilon * eta2 - v

    def tau(v, j):
        return eta2 * hj[j] ** 2 * p1(v) * phi(v, j)

    # bisection starts here
    delta = 0.005
    UB = 999999999
    LB = 0
    while UB - LB > delta:
        v = (float(UB) + LB) / 2
        if Q(v) > 0:
            LB = v
        else:
            UB = v

    x.append(p1(v))
    for j in range(len(M1)):
        x.append(tau(v, j))

    return sum_rate(x), x[0], x[1:]


def lr_method(h):
    # parameters and equations
    o = 100
    p = 3
    u = 0.51
    ki = 10 ** (-2)  # increase the value of ki because the original value is too small
    B = 2 * 10 ** 6
    Vu = 1.1
    N0 = 10 ** (-10)
    wi = np.array([1 if i % 2 == 1 else 1 for i in range(len(h))])  # default weights [1, 1.5, 1, 1.5, 1, 1.5, ...]

    # optimization variables
    tau_i = cvx.Variable(len(h))
    fi = cvx.Variable(len(h))
    ei = cvx.Variable(len(h))

    # optimization objective and constraints
    result = cvx.sum(
        -cvx.multiply(wi, fi) * 10 ** 6 + cvx.multiply(wi, (cvx.kl_div(tau_i, (cvx.multiply(ei, h) + tau_i * N0) / N0) \
                                                            + tau_i - (cvx.multiply(ei,
                                                                                    h) + tau_i * N0) / N0)) * B / Vu / np.log(
            2))
    objective = cvx.Minimize(result)
    constraints = [tau_i >= 0.0, cvx.sum(tau_i) <= 1.0, ei >= 0, fi >= 0,
                   ei + cvx.multiply(ki, fi ** 3) <= u * p * h * (1 - cvx.sum(tau_i))]
    prob = cvx.Problem(objective, constraints)
    rewards = prob.solve(solver=cvx.MOSEK)  # solve the problem by MOSEK

    local_rate = wi * (fi.value) * 10 ** 6
    offloading_rate = wi * B / Vu * tau_i.value * np.log2(1 + ei.value * h / (N0 * tau_i.value))
    mode = []
    for i in range(len(h)):
        if local_rate[i] < offloading_rate[i]:
            mode.append(1)
        else:
            mode.append(0)

    # compute the sum_rate with binary offloading
    E = u * p * h * (1 - np.sum(tau_i.value))
    sum_rate = 0
    for i in range(len(mode)):
        if mode[i] == 1:
            sum_rate += wi[i] * B / Vu * tau_i.value[i] * np.log2(1 + E[i] * h[i] / (N0 * tau_i.value[i]))
        else:
            sum_rate += wi[i] * (E[i] / ki) ** (1 / 3) * 10 ** 6
    # N0 * (math.pow(2, x / B) - 1)
    return sum_rate, mode