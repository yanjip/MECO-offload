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
TB = 1024*GB
PB = 1024*TB

# CPU clock frequency scales
KHZ = 1e3
MHZ = KHZ*1e3
GHZ = MHZ*1e3

# Time scales
slot = 1e-2  # seconds
time_total = 1000  # seconds   total simulation time

# Channels
channel_gain = np.array([-4, -8, -12, -16])  # dB
channel_gain = np.power(10, channel_gain/20)  # multiple
Fading = 1e-3
BW = 10*MHZ  # MHz
N0 = 1e-9  # 单位：W The variance of complex white Gaussian channel noise
TX_POWER = 100  # smallcell=100mWatt
Channel_n = 200  # The number of channels

# MECS
frequency = 20*GHZ  # GHz
f_minportion = 0.1*GHZ  # GHz  # not implemented yet
MEC_CPU_CYCLE = 1500 # cycle per bit; Energy-efficient resource allocation for mobile-edge computation offloading, 2017
MEC_F = 6e9 # cycles per slot.3
MEC_d = 0.2 # VM degrade factor; Reference: A stochastic model to investigate data center performance and qos in iaas cloud computing systems, 2014
DATA = 28597
REQUEST_SIZE =28597
# Reference: Performance and Implications of RAN Caching in LTE Mobile Networks: A Real Traffic Analysis, 2016
# - the total traffic size: 1776GB
# - the total packets: 62104921
# - each packets is 28597 bytes


# Transfer Probability Matrix of UE at any slot
Fk = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],dtype=float)*GHZ  # GHz Unified distributed
Pk = [0,20e-11]               #J/cycle
Rk = [100,500]*KB
Ck = [500, 1500]       #cycles/bit.
power_density = np.array([10, 40])/GHZ**2  # (W/GHz^2) Unified distributed
full_Battery = np.array([1000, 4000])*3.7*3.6
Num_UE=30

# Joule -- mAh*Volt*3.6 Unified distributed
P_send = np.array([3, 8]) # Watt Unified distributed
P_standby = np.array([0.1, 0.2])  # Watt Unified distributed
# Pr(charge_begin) = max(x1 * battery/FullBattery +y1, 0)
# Pr(charge_end) = max(x2 * battery/FullBattery +y2, 0)

# Tasks
Prtask = 0.01 * slot  # Probability of task coming at any slot
data_size = np.array([0.2, 10])*MB  # kB Unified distributed
computation_consumption = np.array([1, 50])*GHZ
