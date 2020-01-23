'''
Jan 2020
A Bayesian phylogenetics inference program
This program is my first practice to understand the concept of phylogenetic reconstruction by using Bayesian
Authors: Nehleh
'''
# import libraries
import sys
import math
from builtins import print

import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt
import seaborn as sns





#=======================================================================================================================
#(A) calculate Maximum Likelihood analytically :

# The first parameter (i.e theta) contains model parameters i.e: in JC69  it just contains edge length v.
# In HKY it contains five parameters the first one is v= edge length, the second one is k = transition-transversion ratio.
# The third last items are pi_A , pi_C and pi_G respectively.
# In GTR: the first one is v= edge length, the third last items are similar to HKY and the four remaining itemes are exchangeability
# parameters
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
        p0 = 0.25 + 0.75 * math.exp(-4 * v / 3)
        p1 = 0.25 - 0.25 * math.exp(-4 * v / 3)
        y = (n-x)* np.log(p0) + x * np.log(p1)
        # negLL = -y
        # print("v = {} , y = {}  , negLL = {} ".format(v, y , negLL))  # for tracing
        return y
        #=====================   for more than two sequences =========================
        # v = theta
        # L = 0.25 * ((0.25 + 0.75 * math.exp(-4 * v[0] / 3)) * (0.25 + 0.75 * math.exp(-4 * v[1] / 3)) * (0.25 - 0.25 * math.exp(-4 * v[2] / 3))
        #  * (0.25 - 0.25 * math.exp(-4 * v[3] / 3)) * (0.25 + 0.75 * math.exp(-4 * v[4] / 3)))
        # negLL = -L
        # return negLL
        #=============================================================================

    elif model == 'HKY':
        ''' '''

    elif model == 'GTR':
        ''' '''
#=======================================================================================================================
# Using scipy.optimize.minimize to find MLE
#=======================================================================================================================
def Max_LL():
    initial_guess = 0.2
    result = spo.minimize(log_likelihood , initial_guess ,method='nelder-mead') #, options={'disp' : True}
    print("Result by using scipy:")
    print("theta = {} , MLE = {}".format(result.x, result.fun))
    return result.x

#=======================================================================================================================
# Using simple loop to find MLE
#=======================================================================================================================
def Max_LL2():
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
    d = np.arange(0.0001, Max_LL() + 0.2, 0.005)
    ll_array = []
    for i in d:
        y_plot = -log_likelihood([i])
        ll_array.append(y_plot)

    plt.plot(d, ll_array, label='fitted model')
    plt.axvline(Max_LL(), color='r', ls='-.')
    plt.title(str(float(Max_LL()))+"  is the maximum likelihood estimate (MLE) of v")
    plt.ylabel('log(L)')
    plt.xlabel('v')
    #plt.legend(loc='lower right')
    plt.show()
#=======================================================================================================================
#The tranistion model defines how to move from current to new
def transition_model(theta):
    propsal = np.random.normal(theta[0], 0.01 ,(1,))
    if ((propsal) < 0):
        return -theta[0]-propsal #reflection
    else:
        return propsal
#=======================================================================================================================
def prior(theta):
        return 1
#=======================================================================================================================
#Defines whether to accept or reject the new sample
def acceptance(current_pos, new_pos):
    if new_pos > current_pos:
        return True
    else:
        accept = np.random.uniform(0,1)
        # Since we did a log likelihood, we need to exponentiate in order to compare to the random number less likely
        # new_pos are less likely to be accepted
        return  (accept < (np.exp(new_pos - current_pos)))
# =======================================================================================================================
def metropolis_hastings(likelihood_computer,prior, transition_model, param_init,iterations,data,acceptance_rule):
    # likelihood_computer(theta,data): returns the likelihood that these parameters generated the data
    # transition_model(theta): a function that draws a sample from a symmetric distribution and returns it
    # param_init: a starting sample
    # iterations: number of accepted to generated
    # data: the data that we wish to model
    # acceptance_rule(current_pos, new_pos): decides whether to accept or reject the new sample


    theta = param_init
    accepted = []
    rejected = []
    final_log_acc = []
    final_log_rej = []
    all_theta = []
    for i in range(iterations):
        new_theta = transition_model(theta)
        all_theta.append(new_theta)
        lik_theta = likelihood_computer(theta,data)
        lik_new_theta = likelihood_computer(new_theta,data)
        if (acceptance(lik_theta+np.log(prior(theta)), lik_new_theta+np.log(prior(new_theta)))):
            theta = new_theta
            accepted.append(new_theta)
            if not(i % 100):
                final_log_acc.append(new_theta)
        else:
            rejected.append(new_theta)
            if not(i % 100):
                final_log_rej.append(new_theta)



    return np.array(final_log_acc) , np.array(final_log_rej) , np.array(accepted) , np.array(rejected) , np.array(all_theta)
# ======================================================================================================================
final_log_acc,final_log_rej,accepted,rejected,all_theta  =  \
    metropolis_hastings(log_likelihood,prior,transition_model,[0.2],1000000,[90,948],acceptance)




print("Accepeted",accepted.shape)

print("Rejected",rejected.shape)



print("Final_Acc",final_log_acc.shape)
print("Final_rej",final_log_rej.shape)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(2,1,1)
ax.plot( rejected[0:500,0], 'rx', label='Rejected',alpha=0.9)
ax.plot( accepted[0:500,0], 'b.', label='Accepted',alpha=0.9)
ax.set_xlabel("Iteration")
ax.set_ylabel("v")
ax.set_title("Figure 2: MCMC sampling for v with Metropolis-Hastings. First 500 samples are shown.")
ax.grid()
ax.legend()

ax2 = fig.add_subplot(2,1,2)
to_show=final_log_acc.shape[0]
ax2.plot( final_log_rej[:to_show,0], 'rx', label='Rejected',alpha=0.5)
ax2.plot( final_log_acc[:to_show,0], 'b.', label='Accepted',alpha=0.5)
ax2.set_xlabel("Iteration")
ax2.set_ylabel("v")
ax2.set_title("Figure 3: MCMC sampling for v with Metropolis-Hastings. All samples are shown.")
ax2.grid()
ax2.legend(loc="best")

fig.tight_layout()


temp=int(0.25*accepted.shape[0])
print("Mean_whole_theta after burn-in: ",np.mean(all_theta[temp:]))
print("SD_whole_theta after burn-in: ",np.std(all_theta[temp:]))

# We consider the initial 25% of the values of v to be "burn-in", so we drop them.
show=int(0.25*final_log_acc.shape[0])

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)
ax.plot(final_log_acc[show:,0])
ax.set_title("Figure 4: Trace for v")
ax.set_ylabel("v")
ax.set_xlabel("Iteration")





fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)
# ax.hist(accepted,bins=300,label="Predicted density")
ax.set_xlabel("$\Theta$")
ax.set_ylabel("Frequency")
ax.set_title("Figure 5: posterior densities for sequence distance $\Theta$ under the JC69 model")
# ax.legend()
ax.grid("off")
ax = sns.distplot(accepted,bins=300,hist= True)

plt.show()


print("This my first change after uploading this file on github")









