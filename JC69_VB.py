import numpy as np
from scipy.stats import gamma
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy.special import digamma

n = 1000
x = 100
alpha = 1
beta = 1


def prior(param):
    d = param
    pr = gamma.pdf(d, alpha)
    return pr


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


def q(param):
    d = param[0]
    gam = param[1]
    gam2 = param[2]
    return gamma.pdf(d, gam, gam2)


def elbo(param):
    d = param[0]
    gam1 = param[1]
    gam2 = param[2]
    eq1 = lambda d: np.log(1 / 4 - 1 / 4 * np.exp(-4 / 3 * d))
    one, err1 = quad(eq1, 0, np.inf)
    eq2 = lambda d: np.log(1 / 4 + 3 / 4 * np.exp(-4 / 3 * d))
    two, err2 = quad(eq2, 0, np.inf)
    three = digamma(gam1) - np.log(gam2)
    total = x * one + (n-x) * two + (alpha - gam1) * three + (gam2 - beta) * gam1/gam2 - gam1 * np.log(gam2) + np.log(gamma(gam1))



    return total


print(elbo(0.1))






# ll_array = []
# d = np.arange(0.01, 0.2, 0.002)
# for i in d:
#     y_plot = posterior_numerical(i)
#     ll_array.append(y_plot)
# plt.plot(d, ll_array, label='true posterior')
# plt.ylabel('density')
# plt.xlabel('d')
# plt.legend(loc='top right')
# plt.show()


