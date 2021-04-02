import numpy as np
from matplotlib import pyplot as plt

# model paramters
T = 1 # final time
m = 3 #	number of assets

# code parameters
n = 1000
h = T/n

# time and BM (each row is an asset)
t = np.r_[1:(n+1)]*h
B = (np.sqrt(h)*np.random.randn(m,n)).cumsum(1)

# stylization parameters; mean and covariance 
mu = np.zeros(m)
sig_a, sig_b = 1., .5
A = sig_a*np.ones((m,m))+(sig_b-sig_a)*np.eye(m)

# solution
S = np.exp((mu-.5*A**2.)*t-A*B)

# wald test
