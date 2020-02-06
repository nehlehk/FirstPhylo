# import libraries
import sys
import math
from builtins import print
from symbol import import_as_name

import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt
import seaborn as sns
import numpy.linalg as la
import pprint


# The first parameter (i.e theta) contains model parameters i.e: in JC69  it just contains edge length v.
#=======================================================================================================================
def log_likelihood(theta, data = None, model = 'JC69'):

    if data == None:
        x = 90
        n = 948
    else:
        x = data[0] # number of mismatch
        n = data[1] # length of alignment

    if model == 'JC69':
        # ===================   for two sequences   =======================
        v = theta[0]
        p0 =  0.25 + 0.75 * math.exp(-4 * v / 3)
        p1 =  0.75 - 0.75 * math.exp(-4 * v / 3)
        y = (n-x)* np.log(p0) + x * np.log(p1)
        # negLL = -y
        # print("v = {} , y = {}  , negLL = {} ".format(v, y , negLL))  # for tracing
        return -y

#=======================================================================================================================
# Using scipy.optimize.minimize to find MLE
#=======================================================================================================================
def max_like():
    initial_guess = 0.2
    result = spo.minimize(log_likelihood , initial_guess ,method='nelder-mead' , options={'disp' : True}) #, options={'disp' : True}
    print("Result by using scipy:")
    print("theta = {} , MLE = {}".format(result.x, -result.fun))
    return result.x

#=======================================================================================================================
# Using simple loop to find MLE
#=======================================================================================================================
def max_like_manual():
    v = np.arange( 0.01, 1, 0.001)
    ll_array = []
    for i in v:
        y_plot = log_likelihood([i])
        ll_array.append(y_plot)

    MLE = np.min(ll_array)
    result = np.where(ll_array == np.amin(ll_array))
    print("Result by using simple loop:")
    print("theta:", float(v[result]), "MLE:", MLE)
    return float(v[result])
#=======================================================================================================================
# Plot Maximum Likelihood Estimation
#=======================================================================================================================
def plot_MLE():
    d = np.arange(0.0001, max_like() + 0.2, 0.005)
    ll_array = []
    for i in d:
        y_plot = -log_likelihood([i])
        ll_array.append(y_plot)

    plt.plot(d, ll_array, label='fitted model')
    plt.axvline(max_like(), color='r', ls='-.')
    plt.title(str(float(max_like()))+"  is the maximum likelihood estimate (MLE) of v")
    plt.ylabel('log(L)')
    plt.xlabel('v')
    #plt.legend(loc='lower right')
    plt.show()

#=======================================================================================================================

# max_like_manual()

plot_MLE()