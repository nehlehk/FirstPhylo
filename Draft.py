# import libraries
import sys
import math
import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt



#=======================================================================================================================
def f(*param):
    x = param[0]
    y = (x - 1.5)**2 + 0.5
    print("x = {} , y = {}". format(x,y)) # for tracing
    return y
#=======================================================================================================================
def test_run():
    initial_guess = 2.0
    min_result = spo.minimize(f , initial_guess , method='nelder-mead', options={'disp' : True} )
    print("Result:")
    print("x = {} , y = {}". format(min_result.x,min_result.fun))

#=======================================================================================================================
def loglik(param, data):
    mu = 0.0
    logLikelihood = 0.0
    resid = 0.0
    print(param[0])
    for i in range(data.shape[0]):
        mu = param[0] * (data[i][0] ** param[1])
        resid = float(data[i][1] - mu) / param[2]
        logLikelihood += -0.5 * resid * resid

    return -(-data.shape[0] * np.log(param[2]) + logLikelihood)
#=======================================================================================================================
def Max_LL2():
    d = np.arange( - 5, 5, 0.5)
    ll_array = []
    for i in d:
        y_plot = f(i)
        ll_array.append(y_plot)

    MLE =np.min(ll_array)
    result = np.where(ll_array == np.amin(ll_array))
    print("Min:",float(d[result]), "Min value:",MLE)
#=======================================================================================================================

#test_run()
# loglik(3,7,[3,4,5])
#Max_LL2()
generate_samples()