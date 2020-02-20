import numpy as np
from scipy.stats import gamma
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy.special import digamma
from scipy.special import loggamma
from scipy.optimize import Bounds
import scipy.optimize as spo

n = 1000
x = 100
alpha = 1
beta = 1



def prior(param):
    return gamma(a=alpha, scale=beta, loc=0).pdf(param)



def p(param):
    d = param
    p0 = 1 / 4 + 3 / 4 * np.exp(-4 / 3 * d)
    p1 = 1 / 4 - 1 / 4 * np.exp(-4 / 3 * d)
    res = p0 ** (n - x) * p1 ** x
    return  res


def posterior_numerical(param):
    numer = prior(param) * p(param)
    temp_denom = lambda d: prior(d) * p(d)
    denom , err = quad(temp_denom, 0, np.inf)
    post = numer / denom
    return post



def elbo(param):
    d = param[0]
    gam1 = param[1]
    gam2 = param[2]
    eq1 = lambda d: np.log(1 / 4 - 1 / 4 * np.exp(-4 / 3 * d))
    one, err1 = quad(eq1, 0, np.inf)
    eq2 = lambda d: np.log(1 / 4 + 3 / 4 * np.exp(-4 / 3 * d))
    two, err2 = quad(eq2, 0, np.inf)
    three = digamma(gam1) - np.log(gam2)
    total = (x * one)  \
            + ((n-x) * two)  \
            + (alpha - gam1) * three \
            + (gam2 - beta) * (gam1/gam2) \
            - gam1 * np.log(gam2) \
            + loggamma(gam1)
    return total



def max_param():
    initial_guess = [0.2,1,1]
    bounds =  Bounds([0 , 0 , 0 ], [1, 2, 2])
    result = spo.minimize(elbo , initial_guess , method='trust-constr' , bounds = bounds, options={'verbose': 1} )
    print("Result by using scipy:")
    print("edge_length = {} , gamma1 ={} , gamma2 = {}".format(result.x[0],result.x[1],result.x[2]))
    return result



max_param()

ll_array = []
d = np.arange(0.01, 0.2, 0.002)
for i in d:
    y_plot = posterior_numerical(i)
    ll_array.append(y_plot)
plt.plot(d, ll_array, label='true posterior')
plt.ylabel('density')
plt.xlabel('d')
plt.legend(loc='upper right')
plt.show()


