# import libraries
import numpy as np
import scipy.optimize as spo
import pprint

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
def count_v1(dna,base):
    # i = 0
    # for c in dna:
    #     if c== base:
    #         i += 1
    # return i
     return dna.count(base)
#====================================================================================================================
def freq_lists(dna_list):
    n = len(dna_list[0])
    A = [0] * n
    C = [0] * n
    T = [0] * n
    G = [0] * n
    for dna in dna_list:
        for index , base in enumerate(dna):
            if base == 'A':
                A[index] += 1
            elif base == 'C':
                C[index] += 1
            elif base == 'G':
                G[index] += 1
            elif base == 'T':
                T[index] += 1


    return A,C,G,T
#===================================================================================================================
def frequency_matrix_computation(dna_list):
    n = len(dna_list[0])
    frequency_matrix = {base: np.zeros(n , dtype=np.int)  for base in 'ACGT'}

    for dna in dna_list:
        # dna = np.array(dna , dtype = 'c')
        for index , base in enumerate(dna):
            frequency_matrix[base][index] +=1

    return  frequency_matrix
#===================================================================================================================
def get_base_count(dna):
    counts = {'A':0 , 'C':0 , 'G':0 , 'T':0}
    for base in dna:
        counts[base] +=1

    return counts
#===================================================================================================================
def get_base_frequencies(dna):
    # base_count = {base: dna.count(base)   for base in 'ACGT'}
    base_freq =  {base: dna.count(base)/float(len(dna))
             for base in 'ACGT'}

    return  base_freq
#===================================================================================================================
def transition_matrix_computation(seq1,seq2):
    n = len(seq1)
    tran_matrix = np.zeros((4,4) , dtype= int )
    for i in range(n):
      if seq1[i] == seq2[i] == 'A':
          tran_matrix[0][0] +=1
      elif seq1[i] == seq2[i] == 'C':
          tran_matrix[1][1] +=1
      elif seq1[i] == seq2[i] == 'G':
          tran_matrix[2][2] += 1
      elif seq1[i] == seq2[i] == 'T':
          tran_matrix[3][3] += 1
      elif seq1[i] == 'A' and seq2[i] == 'C':
          tran_matrix[0][1] += 1
      elif seq1[i] == 'A' and seq2[i] == 'G':
          tran_matrix[0][2] += 1
      elif seq1[i] == 'A' and seq2[i] == 'T':
          tran_matrix[0][3] += 1
      elif seq1[i] == 'C' and seq2[i] == 'A':
          tran_matrix[1][0] += 1
      elif seq1[i] == 'C' and seq2[i] == 'G':
          tran_matrix[1][2] += 1
      elif seq1[i] == 'C' and seq2[i] == 'T':
          tran_matrix[1][3] += 1
      elif seq1[i] == 'G' and seq2[i] == 'A':
          tran_matrix[2][0] += 1
      elif seq1[i] == 'G' and seq2[i] == 'C':
          tran_matrix[2][1] += 1
      elif seq1[i] == 'G' and seq2[i] == 'T':
          tran_matrix[2][3] += 1
      elif seq1[i] == 'T' and seq2[i] == 'A':
          tran_matrix[3][0] += 1
      elif seq1[i] == 'T' and seq2[i] == 'C':
          tran_matrix[3][1] += 1
      elif seq1[i] == 'T' and seq2[i] == 'G':
          tran_matrix[3][2] += 1

    return tran_matrix

#===================================================================================================================
def transition_matrix_computation2(dna_list):
    i = 0
    n = len(dna_list)
    tran_matrix = {base :{base2 :0 for base2 in 'ACGT'} for base in 'ACGT'}

    # while (i < n-1 ):
    dna1 = dna_list[i]
    dna2 = dna_list[i+1]
    for index, base1 in enumerate(dna1):
        base2 = dna2[index]
        tran_matrix[base1][base2] += 1

        # i += 1

    return tran_matrix
#===================================================================================================================
def avg_base_frequencies(alignment):
    n = len(alignment)
    # sum_bf = [0,0,0,0]
    for align in alignment:
       sum_bf =  get_base_frequencies(align)

    return sum_bf
#===================================================================================================================
fr = avg_base_frequencies(['ACCT','ACGT'])
print(fr['A']/2)
# tran_matrix2 = {base :{base2 :0 for base2 in 'ACGT'} for base in 'ACGT'}
#
# tran_matrix2['A']['A'] = tran_matrix2['A']['A'] + 1
# pprint.pprint(tran_matrix2)

# t = transition_matrix_computation2\
#     (['GAAGTCCTTGAGAAATAAACTGCACACACTGG','GGACTCCTTGAGAAATAAACTGCACACACTGG','GGACTCCTTGAGAAATAAACTGCACACACTGG'])
# pprint.pprint(t)


# print(tran_matrix['A']['A'])

# pprint.pprint(frequency_matrix_computation(['GAAGTCCTTGAGAAATAAACTGCACACACTGG','GGACTCCTTGAGAAATAAACTGCACACACTGG']))
# t = transition_matrix_computation('GAAGTCCTTGAGAAATAAACTGCACACACTGG','GGACTCCTTGAGAAATAAACTGCACACACTGG')
# pprint.pprint(t)
# tran_matrix = np.zeros((4,4) , dtype= int )
#
# print(tran_matrix[3][3])


# print(np.random.choice(list('ACGT')))

# pprint.pprint(get_base_frequencies('ACCCCTGG'))
# c , f = get_base_frequencies('ACCCCTGG')
#
# print(c['C'])
# print(f['C'])
# print(get_base_count('ACCCCTGG'))

# frequency_matrix = frequency_matrix_computation(['ACCCT','ACGTC'])

# print(frequency_matrix['A'][0])
#
# pprint.pprint(frequency_matrix)

# A, C, G, T = freq_lists(['GAAG','GGAC'])
# print(A)
# print(C)
# print(G)
# print(T)



# result = count_v1('GAAGTCCTTGAGAAATAAACTGCACACACTGG','G')
# print(result)

#test_run()
# loglik(3,7,[3,4,5])
#Max_LL2()
