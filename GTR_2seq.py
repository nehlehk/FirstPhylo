#======================================================================================================================

#======================================================================================================================
# import libraries

from builtins import print
import numpy as np
import scipy.optimize as spo
from scipy import optimize
import numpy.linalg as la
import pprint
import dendropy
from scipy.optimize import Bounds

#===================================================================================================================
def get_base_frequencies(dna):
    base_freq =  {base: dna.count(base)/float(len(dna)) for base in 'ACGT'}

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
def give_index(c):
    if c == "A":
        return 0
    elif c == "C":
        return 1
    elif c == "G":
        return 2
    elif c == "T":
        return 3
#=======================================================================================================================
def cons_rate(x):
    return x[1] + x[2] + x[3] + x[4] + x[5]  - 1
#=======================================================================================================================
def GTR(params):
    n = alignment.sequence_size
    mu = 0

    freq = np.zeros((4, 4))
    q = np.zeros((4, 4))
    sqrtPi =  np.zeros((4, 4))
    sqrtPiInv =  np.zeros((4, 4))
    exchang = np.zeros((4, 4))
    s = np.zeros((4, 4))
    fun = np.zeros(n)
    a,b,c,d,e  = params[1:]
    f = 1

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

    p = np.dot(left, np.dot(np.diag(np.exp(eigval * params[0])), right))

    k = 0
    dna1 = alignment[k]
    dna2 = alignment[k + 1]
    for index, base1 in enumerate(dna1):
       base2 = dna2[index]
       i = give_index(str(base1))
       j = give_index(str(base2))
       fun[index] = pi[i] * p[i][j]

    ll = np.sum(np.log(fun))

    pprint.pprint(p)

    return -ll
#===================================================================================================================
def constraint(rate):
    return sum(rate) - 1
#=======================================================================================================================
def max_likeGTR():
    initial_guess = [0.1, 0,0.5,0.5,0.5,0.5]
    cons = {"type": "eq", "fun": cons_rate}
    bounds =  Bounds([0 , 0 , 0 , 0 , 0 , 0 ], [1000, 1000, 1000, 1000, 1000, 1000])
    result = spo.minimize(GTR , initial_guess , method='COBYLA' , options={'verbose': 1}, bounds = bounds ) #,constraints= cons
    print("Result by using scipy:")
    print("theta = {} , MLE = {}".format(result.x, -result.fun))
    return result
#=======================================================================================================================



# a= b = c = d = e =  f = 1

pi = [0.25]*4


alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/0_Research/PhD/Data/test.fasta"), schema="fasta")

max_likeGTR()


