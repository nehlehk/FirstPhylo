import numpy as np
from hmmlearn import hmm
import matplotlib.pyplot as plt
import pandas as pd



df1 = pd.read_csv("/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/faketree/RAxML_perSiteLLs.likelihood_GTR", sep='\s+', header=None)
ll = df1.to_numpy()
data = np.array(ll)
X = data.reshape((-1,1))
X = X[1:100000]

mean = np.mean(X)
std = np.std(X)
# print(mean)
# print(std)

# a = float(np.random.normal(mean,std,1))
# b = float(np.random.normal(mean,std,1))

a= -6.970009569914
b= mean
astd= std
bstd= std
print(a)
print(b)




model = hmm.GaussianHMM(n_components=2, covariance_type="full" ,algorithm='viterbi' )
model.startprob_ = np.array([0.88, 0.12])
# model.startprob_ = np.array([0.9, 0.1])
model.transmat_ = np.array([[0.9999, 0.0001] , [0.0001, 0.9999]])
model.means_ = np.array([[a, astd], [b, bstd]])
model.covars_ = np.tile(np.identity(2), (2, 1, 1))


posterior = model.predict_proba(X)

hiddenStates = model.predict(X)
print(hiddenStates)

score = model.score(X)

fig = plt.figure(figsize=(15,8))
ax = fig.add_subplot(2,1,1)
ax.set_title("Hidden Markov Models - ClonalFrame and Recombination -- log probability of the most likely state is  " + str (score))
ax.plot(hiddenStates)
ax.set_ylabel("Clonal - NonClonal State")



ax2 = fig.add_subplot(2,1,2)
ax2.plot(posterior)
ax2.set_ylabel("posterior probability for each state")
plt.show()












# newmodel =  hmm.GaussianHMM(n_components=2).fit(X)
#
# print("Transitin:")
# print(newmodel.transmat_)
# print("emission:")
# print(newmodel.means_)
#
# posterior = newmodel.predict_proba(X)
# hiddenStates = newmodel.predict(X)
# score = newmodel.score(X)
#
# fig = plt.figure(figsize=(15,8))
# ax = fig.add_subplot(2,1,1)
# ax.set_title("Hidden Markov Models - ClonalFrame and Recombination -- log probability of the most likely state is  " + str (score))
# ax.plot(hiddenStates)
# ax.set_ylabel("Clonal - NonClonal State")
#
#
# ax2 = fig.add_subplot(2,1,2)
# ax2.plot(posterior)
# ax2.set_ylabel("posterior probability for each state")
# plt.show()
#
