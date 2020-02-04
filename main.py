#======================================================================================================================
'''
Jan 2020
Bayesian phylogenetics inference code
This code is my first practice to understand the concept of phylogenetic reconstruction by using Bayesian.
Authors: Nehleh Kargarfard
'''
#======================================================================================================================
# import libraries
import sys
import math
from builtins import print
import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt
import seaborn as sns
import numpy.linalg as la
import pprint
# from Bio import SeqIO


#=======================================================================================================================
def transition_probability_matrix_GTR(alignment,br_len ,a = 1 , b = 1 , c = 1 , d = 1 , e = 1 , f = 1):
    n = len(alignment[0])
    mu = 0

    freq = np.zeros((4, 4))
    q = np.zeros((4, 4))
    sqrtPi =  np.zeros((4, 4))
    sqrtPiInv =  np.zeros((4, 4))
    exchang = np.zeros((4, 4))
    s = np.zeros((4, 4))

    # pi = avg_base_frequencies(alignment)
    pi = [0.25,0.25,0.25,0.25]
    freq = np.diag(pi)
    sqrtPi = np.diag(np.sqrt(pi))
    sqrtPiInv = np.diag(1.0/np.sqrt(pi))

    mu = 1 / (2 * (a * pi[0] * pi[1]) + (b * pi[0] * pi[2]) + (c * pi[0] * pi[3]) + (d * pi[1] * pi[2]) + (
                e * pi[1] * pi[3]) + (pi[2] * pi[3]))

    exchang[0][1] = exchang[1][0] = a
    exchang[0][2] = exchang[2][0] = b
    exchang[0][3] = exchang[3][0] = c
    exchang[1][2] = exchang[2][1] = d
    exchang[1][3] = exchang[3][1] = e
    exchang[2][3] = exchang[3][2] = f

    q = np.multiply(np.dot(exchang,freq), mu)


    for i in range(4):
        q[i][i] = -sum(q[i][0:4])


    s = np.dot(sqrtPi,np.dot(q,sqrtPiInv))

    eigval, eigvec = la.eig(s)
    eigvec_inv = la.inv(eigvec)

    left = np.dot(sqrtPi,eigvec)
    right = np.dot(eigvec_inv,sqrtPiInv)

    p = np.dot(left, np.dot(np.diag(np.exp(eigval * br_len)), right))

    # i = 0
    # dna1 = alignment[i]
    # dna2 = alignment[i + 1]
    # for index, base1 in enumerate(dna1):
    #    base2 = dna2[index]
    #    f[index] =  p[give_index(base1)][give_index(base1)] * p[give_index(base1)][give_index(base1)]
    #
    #
    # result = np.log(np.sum(f))
    #
    # print(result)

    return p
#===================================================================================================================
def give_index(c):
    if c == "A":
        return 0
    elif c == "C":
        return 1
    elif c == "G":
        return 2
    elif c == "T":
        return 3
#===================================================================================================================
def ratio_matrix_computation(dna_list):
    i = 0
    seq_count = len(dna_list)
    seq_length = len(dna_list[0])
    ratio_matrix = {base :{base2 :0 for base2 in 'ACGT'} for base in 'ACGT'}

    # while (i < seq_count -1 ):
    dna1 = dna_list[i]
    dna2 = dna_list[i+1]
    for index, base1 in enumerate(dna1):
        base2 = dna2[index]
        ratio_matrix[base1][base2] += 1  #give count

    # ============================================================
    #dirty code :(
    for base1 in 'ACGT':
        for base2 in 'ACGT':
            ratio_matrix[base1][base2] = ratio_matrix[base1][base2] / seq_length # give ratio

        # i += 1

    return ratio_matrix
#===================================================================================================================
def get_base_frequencies(dna):
    base_freq =  {base: dna.count(base)/float(len(dna))
             for base in 'ACGT'}

    return  base_freq
#===================================================================================================================
def avg_base_frequencies(alignment):
    n = len(alignment)
    a = c = g = t = 0
    for align in alignment:
       sum_bf =  get_base_frequencies(align)
       a += sum_bf['A']
       c += sum_bf['C']
       g += sum_bf['G']
       t += sum_bf['T']

    return [a/float(n) , c/float(n) , g/float(n) , t/float(n)]
#===================================================================================================================
def  transition_probability_matrix_JC(t,beta):
    q = np.zeros((4, 4))
    p = np.zeros((4, 4))
    q[:][:] = beta
    q[0][0] = q[1][1] = q[2][2] = q[3][3] = -3*beta
    eigval, eigvec = la.eig(q)
    eigvec_inv = la.inv(eigvec)
    p = np.dot(eigvec,np.dot(np.diag(np.exp(eigval * t)), eigvec_inv))
    return p
#=======================================================================================================================
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
        v = theta
        p0 =  0.25 + 0.75 * math.exp(-4 * v / 3)
        p1 =  0.75 - 0.75 * math.exp(-4 * v / 3)
        y = (n-x)* np.log(p0) + x * np.log(p1)
        # negLL = -y
        # print("v = {} , y = {}  , negLL = {} ".format(v, y , negLL))  # for tracing
        return y
        #=====================   for more than two sequences =========================

        #=============================================================================

    elif model == 'HKY':
        ''' '''

    elif model == 'GTR':
        ''' '''
#=======================================================================================================================
# Using scipy.optimize.minimize to find MLE
#=======================================================================================================================
def max_like():
    initial_guess = 0.2
    result = spo.minimize(log_likelihood , initial_guess ,method='nelder-mead') #, options={'disp' : True}
    print("Result by using scipy:")
    print("theta = {} , MLE = {}".format(result.x, result.fun))
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
#The tranistion model defines how to move from current to new
def transition_model(theta,type,width):
    if type == 'sw_norm':
        proposal = np.random.normal(theta, width)
        if (proposal < 0):
            proposal =  -theta - proposal  # reflection in normal dist --- the excess is reflected back into the interval
    elif type == 'sw_uni':
        proposal = np.random.uniform(theta - width/2 , theta + width/2)
        if (proposal < 0):
            proposal = - proposal  # reflection in uniform dist
    elif type == 'pscale':
        c = np.exp(width * (np.random.uniform(0,1) -1)/2 )
        proposal = theta * c


    # print(proposal)
    return proposal
#=======================================================================================================================
#Define prior
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


    theta = param_init[0]
    proposal_type = param_init[1]
    proposal_width = param_init[2]
    accepted = []
    rejected = []
    final_log_acc = []
    final_log_rej = []
    all_theta = []
    for i in range(iterations):
        new_theta = transition_model(theta,proposal_type,proposal_width)
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
        all_theta.append(theta)


    return np.array(final_log_acc) , np.array(final_log_rej) , np.array(accepted) , np.array(rejected) , np.array(all_theta)
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
v = np.arange( 0 , 2, 0.01)
y1 = []
y2 = []
d = []
for i in v:
    # p = transition_probability_matrix_JC(i,i)
    p = transition_probability_matrix_GTR(['AAAAAGGCAA', 'GGGCTCTTAA'],i)
    y1.append(p[0][0])
    y2.append(p[0][1])
    d.append(i * i)



pprint.pprint(p)
plt.plot(d, y1, label='P(AA)')
plt.plot(d, y2, label='P(AC)')
plt.ylabel('probablity')
plt.xlabel('v')
plt.axhline(p[0][0], color='r', ls='-.')
plt.legend(loc='lower right')
plt.show()




# iterations = 1000000
#
# # transition_model(0.5,'pscale',0.1)
#
# final_log_acc,final_log_rej,accepted,rejected,all_theta  =  \
#     metropolis_hastings(log_likelihood,prior,transition_model,[0.5,'sw_uni',0.15],iterations,[90,948],acceptance)
#
#
#
# print("P_jump = ", round(accepted.shape[0]/iterations * 100 , 2) , '%' )
#
# # print("Rejected rate",int(rejected.shape)/iterations)
#
#
#
# print("Final_Acc",final_log_acc.shape)
# print("Final_rej",final_log_rej.shape)
#
# fig = plt.figure(figsize=(10,10))
# ax = fig.add_subplot(2,1,1)
# ax.plot( rejected[:500], 'rx', label='Rejected',alpha=0.9)
# ax.plot( accepted[:500], 'b.', label='Accepted',alpha=0.9)
# ax.set_xlabel("Iteration")
# ax.set_ylabel("v")
# ax.set_title("Figure 2: MCMC sampling for v with Metropolis-Hastings. First 500 samples are shown.")
# ax.grid()
# ax.legend()
#
# ax2 = fig.add_subplot(2,1,2)
# to_show=final_log_acc.shape[0]
# ax2.plot( final_log_rej[:to_show], 'rx', label='Rejected',alpha=0.5)
# ax2.plot( final_log_acc[:to_show], 'b.', label='Accepted',alpha=0.5)
# ax2.set_xlabel("Iteration")
# ax2.set_ylabel("v")
# ax2.set_title("Figure 3: MCMC sampling for v with Metropolis-Hastings. All samples are shown.")
# ax2.grid()
# ax2.legend(loc="best")
#
# fig.tight_layout()
#
#
# temp=int(0.25*all_theta.shape[0])
# print("Mean_whole_theta after burn-in: ",np.mean(all_theta[temp:]))
# print("SD_whole_theta after burn-in: ",np.std(all_theta[temp:]))
#
# # consider the initial 25% of the values of v to be "burn-in", so we drop them.
# show=int(0.25*final_log_acc.shape[0])
#
# fig = plt.figure(figsize=(10,10))
# ax = fig.add_subplot(1,1,1)
# ax.plot(final_log_acc[show:])
# ax.set_title("Figure 4: Trace for v")
# ax.set_ylabel("v")
# ax.set_xlabel("Iteration")
#
#
#
#
#
# fig = plt.figure(figsize=(10,10))
# ax = fig.add_subplot(1,1,1)
# # ax.hist(accepted,bins=300,label="Predicted density")
# ax.set_xlabel("$\Theta$")
# ax.set_ylabel("Frequency")
# ax.set_title("Figure 5: posterior densities for sequence distance $\Theta$ under the JC69 model")
# # ax.legend()
# ax.grid("off")
# ax = sns.distplot(accepted,bins=300,hist= True)
#
# plt.show()

