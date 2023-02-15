'''
@Author  ：Yan JP
@Created on Date：2023/2/4 14:59 
'''
from numpy import random
import numpy as np
x = random.rayleigh(scale=1, size=(3,3))
print('标准差为2，大小为3×3的瑞利分布：\n', x)

# 可以取整
R = np.trunc((1500 + 300*(np.random.uniform(0,1,5))))
R1=R.view()
print(R1)
for i in R1:
    print(i)
print(type(R1))
a=np.ones(5)*10
print(a)

from scipy.special import lambertw
print(lambertw(1))

import cvxpy
import math
print(math.sqrt(np.random.exponential(1e-6)))

lk=[0,0,0]
Rk=[1,2,4]
lk[1:3]=Rk[1:3]
print(lk)

from sympy import *
