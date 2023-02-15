'''
@Author  ：Yan JP
@Created on Date：2023/2/3 22:42 
'''
import math
import numpy as np

# Data size scales
BYTE = 8
KB = 1024*BYTE
MB = 1024*KB
GB = 1024*MB
# TB = 1024*GB
# PB = 1024*TB

# CPU clock frequency scales
KHZ = 1e3
MHZ = KHZ*1e3
GHZ = MHZ*1e3

# Time scales
T_slot =1.05# seconds  0.5-0.8下降曲线（0.8的时候特别小）      1.15

#  UE at any slot
Fk = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],dtype=float)*GHZ  # GHz Unified distributed
Pk = [0,20e-11]               #J/cycle
# Rk = [100*KB,500*KB]
Rk = [100,500]

Ck = [500, 1500]       #cycles/bit.
UE_n=25


# Channels
B = 10*MHZ  # MHz
N0 = 1e-9  # 单位：W The variance of complex white Gaussian channel noise

# MECS
F_MEC =61e9 # cycles per slot.3 33970851290.0  57-69      69-28-2.3







import  numpy as np
import math


def fx(x):
    return N0*(math.pow(2,x/B)-1)
