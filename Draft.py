# import libraries
import numpy
import scipy.optimize as spo
import pprint
import numpy.linalg as la
import dendropy
import math



alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/0_Research/PhD/Data/LL_vector/JC69_100.fasta"), schema="fasta")

for i in range(alignment.sequence_size):
    s = ''
    for j in range(len(alignment)):
        s += str(alignment[j][i])
        # for idx, name in enumerate(alignment[i][:]):
        #     print(idx,name)
    print(s)


# for d in range(alignment.sequence_size):
#     print(alignment[][d])
    # for idx, name in enumerate(alignment[d][:]):
    #     print(idx,name)




# print(ll)

# print(partial)



# def give_index(c):
#     if c == "A":
#         return 0
#     elif c == "C":
#         return 1
#     elif c == "G":
#         return 2
#     elif c == "T":
#         return 3
# # =======================================================================================================
# column = 'TCACC'
#
# s = 5 #tips number
# br_length = 0.1
# p = np.zeros((4, 4))
#
# mytree = dendropy.Tree.get_from_string('(((1,2),3),(4,5));', 'newick')
# partial = np.zeros((4,2*s -2))
#
# for i in range(4):
#     for j in range(4):
#         p[i][j] = 0.25 - 0.25* math.exp(-4*br_length/3)
#     p[i][i] = 0.25 + 0.75*math.exp(-4*br_length/3)
# pprint.pprint(p)
#
# post_order =[]
# for  index, lf in enumerate(mytree.leaf_nodes()):
#     i = give_index(column[index])
#     partial[i][index] = 1
#
#
#
# for index,node in enumerate(mytree.postorder_node_iter()):
#     post_order.append(node.taxon)
#     # if node.is_leaf():
#     #     i = give_index(column[index])
#     #     partial[i][index] = 1
#     # else:
#     if not node.is_leaf():
#         children = node.child_nodes()
#         leftchild =  post_order.index(children[0].taxon)
#         rightchild = post_order.index(children[1].taxon)
#         for j in range(4):
#             left = 0
#             right = 0
#             for k in range(4):
#                 left += p[j][k] * partial[k][leftchild]
#                 right += p[j][k] * partial[k][rightchild]
#             partial[j][index] = left * right
#
#
#
# print(partial)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# # alignment = dendropy.DnaCharacterMatrix.get\
# #     (file=open("/home/nehleh/0_Research/PhD/Data/test.fasta"), schema="fasta" , )
# # print(alignment.sequence_size)
# # print(len(alignment))
# # for seq in alignment.values():
# #    print(seq)
# #=======================================================================================================================
# def f(*param):
#     x = param[0]
#     y = (x - 1.5)**2 + 0.5
#     print("x = {} , y = {}". format(x,y)) # for tracing
#     return y
# #=======================================================================================================================
# def test_run():
#     initial_guess = 2.0
#     min_result = spo.minimize(f , initial_guess , method='nelder-mead', options={'disp' : True} )
#     print("Result:")
#     print("x = {} , y = {}". format(min_result.x,min_result.fun))
#
# #=======================================================================================================================
# def loglik(param, data):
#     mu = 0.0
#     logLikelihood = 0.0
#     resid = 0.0
#     print(param[0])
#     for i in range(data.shape[0]):
#         mu = param[0] * (data[i][0] ** param[1])
#         resid = float(data[i][1] - mu) / param[2]
#         logLikelihood += -0.5 * resid * resid
#
#     return -(-data.shape[0] * np.log(param[2]) + logLikelihood)
# #=======================================================================================================================
# def Max_LL2():
#     d = np.arange( - 5, 5, 0.5)
#     ll_array = []
#     for i in d:
#         y_plot = f(i)
#         ll_array.append(y_plot)
#
#     MLE =np.min(ll_array)
#     result = np.where(ll_array == np.amin(ll_array))
#     print("Min:",float(d[result]), "Min value:",MLE)
# #=======================================================================================================================
# def count_v1(dna,base):
#     # i = 0
#     # for c in dna:
#     #     if c== base:
#     #         i += 1
#     # return i
#      return dna.count(base)
# #====================================================================================================================
# def freq_lists(dna_list):
#     n = len(dna_list[0])
#     A = [0] * n
#     C = [0] * n
#     T = [0] * n
#     G = [0] * n
#     for dna in dna_list:
#         for index , base in enumerate(dna):
#             if base == 'A':
#                 A[index] += 1
#             elif base == 'C':
#                 C[index] += 1
#             elif base == 'G':
#                 G[index] += 1
#             elif base == 'T':
#                 T[index] += 1
#
#
#     return A,C,G,T
# #===================================================================================================================
# def frequency_matrix_computation(dna_list):
#     n = len(dna_list[0])
#     frequency_matrix = {base: np.zeros(n , dtype=np.int)  for base in 'ACGT'}
#
#     for dna in dna_list:
#         # dna = np.array(dna , dtype = 'c')
#         for index , base in enumerate(dna):
#             frequency_matrix[base][index] +=1
#
#     return  frequency_matrix
# #===================================================================================================================
# def get_base_count(dna):
#     counts = {'A':0 , 'C':0 , 'G':0 , 'T':0}
#     for base in dna:
#         counts[base] +=1
#
#     return counts
# #===================================================================================================================
# def get_base_frequencies(dna):
#     base_freq =  {base: dna.count(base)/float(len(dna))
#              for base in 'ACGT'}
#
#     return  base_freq
# #===================================================================================================================
# def transition_matrix_computation(seq1,seq2):
#     n = len(seq1)
#     tran_matrix = np.zeros((4,4) , dtype= int )
#     for i in range(n):
#       if seq1[i] == seq2[i] == 'A':
#           tran_matrix[0][0] +=1
#       elif seq1[i] == seq2[i] == 'C':
#           tran_matrix[1][1] +=1
#       elif seq1[i] == seq2[i] == 'G':
#           tran_matrix[2][2] += 1
#       elif seq1[i] == seq2[i] == 'T':
#           tran_matrix[3][3] += 1
#       elif seq1[i] == 'A' and seq2[i] == 'C':
#           tran_matrix[0][1] += 1
#       elif seq1[i] == 'A' and seq2[i] == 'G':
#           tran_matrix[0][2] += 1
#       elif seq1[i] == 'A' and seq2[i] == 'T':
#           tran_matrix[0][3] += 1
#       elif seq1[i] == 'C' and seq2[i] == 'A':
#           tran_matrix[1][0] += 1
#       elif seq1[i] == 'C' and seq2[i] == 'G':
#           tran_matrix[1][2] += 1
#       elif seq1[i] == 'C' and seq2[i] == 'T':
#           tran_matrix[1][3] += 1
#       elif seq1[i] == 'G' and seq2[i] == 'A':
#           tran_matrix[2][0] += 1
#       elif seq1[i] == 'G' and seq2[i] == 'C':
#           tran_matrix[2][1] += 1
#       elif seq1[i] == 'G' and seq2[i] == 'T':
#           tran_matrix[2][3] += 1
#       elif seq1[i] == 'T' and seq2[i] == 'A':
#           tran_matrix[3][0] += 1
#       elif seq1[i] == 'T' and seq2[i] == 'C':
#           tran_matrix[3][1] += 1
#       elif seq1[i] == 'T' and seq2[i] == 'G':
#           tran_matrix[3][2] += 1
#
#     return tran_matrix
#
# #===================================================================================================================
# def ratio_matrix_computation(dna_list):
#     i = 0
#     seq_count = len(dna_list)
#     seq_length = len(dna_list[0])
#     ratio_matrix = {base :{base2 :0 for base2 in 'ACGT'} for base in 'ACGT'}
#
#     # while (i < seq_count -1 ):
#     dna1 = dna_list[i]
#     dna2 = dna_list[i+1]
#     for index, base1 in enumerate(dna1):
#         base2 = dna2[index]
#         ratio_matrix[base1][base2] += 1  #give count
#
#     # ============================================================
#     #dirty code :(
#     for base1 in 'ACGT':
#         for base2 in 'ACGT':
#             ratio_matrix[base1][base2] = ratio_matrix[base1][base2] / seq_length # give ratio
#
#         # i += 1
#
#     return ratio_matrix
# #===================================================================================================================
# def avg_base_frequencies(alignment):
#     n = len(alignment)
#     a = c = g = t = 0
#     for align in alignment:
#        sum_bf =  get_base_frequencies(align)
#        a += sum_bf['A']
#        c += sum_bf['C']
#        g += sum_bf['G']
#        t += sum_bf['T']
#
#     return [a/float(n) , c/float(n) , g/float(n) , t/float(n)]
# #===================================================================================================================
# def transition_probability_matrix_GTR(alignment,t):
#     n= alignment.sequence_size
#     # n = len(alignment[0])
#     mu = 0
#     a = b = c = d = e = f = 1
#
#     freq = np.zeros((4, 4))
#     q = np.zeros((4, 4))
#     sqrtPi =  np.zeros((4, 4))
#     sqrtPiInv =  np.zeros((4, 4))
#     exchang = np.zeros((4, 4))
#     s = np.zeros((4, 4))
#
#     # ratio_matrix =ratio_matrix_computation(alignment)
#
#     # pi = avg_base_frequencies(alignment)
#     pi = [0.25,0.25,0.25,0.25]
#     freq = np.diag(pi)
#     sqrtPi = np.diag(np.sqrt(pi))
#     sqrtPiInv = np.diag(1.0/np.sqrt(pi))
#
#     mu = 1 / (2 * (a * pi[0] * pi[1]) + (b * pi[0] * pi[2]) + (c * pi[0] * pi[3]) + (d * pi[1] * pi[2]) + (
#                 e * pi[1] * pi[3]) + (pi[2] * pi[3]))
#
#     exchang[0][1] = exchang[1][0] = a
#     exchang[0][2] = exchang[2][0] = b
#     exchang[0][3] = exchang[3][0] = c
#     exchang[1][2] = exchang[2][1] = d
#     exchang[1][3] = exchang[3][1] = e
#     exchang[2][3] = exchang[3][2] = f
#
#     q = np.multiply(np.dot(exchang,freq), mu)
#
#
#     for i in range(4):
#         q[i][i] = -sum(q[i][0:4])
#
#
#     s = np.dot(sqrtPi,np.dot(q,sqrtPiInv))
#
#     eigval, eigvec = la.eig(s)
#     eigvec_inv = la.inv(eigvec)
#
#     left = np.dot(sqrtPi,eigvec)
#     right = np.dot(eigvec_inv,sqrtPiInv)
#
#     p = np.dot(left, np.dot(np.diag(np.exp(eigval * t)), right))
#
#     # pprint.pprint(p)
#     #
#     return p
# #===================================================================================================================
# def  transition_probability_matrix_JC(t,beta):
#     q = np.zeros((4, 4))
#     p = np.zeros((4, 4))
#     q[:][:] = beta
#     q[0][0] = q[1][1] = q[2][2] = q[3][3] = -3*beta
#     # pprint.pprint(q)
#     eigval, eigvec = la.eig(q)
#     eigvec_inv = la.inv(eigvec)
#     p = np.dot(eigvec,np.dot(np.diag(np.exp(eigval * t)), eigvec_inv))
#     return p
# #===================================================================================================================
# def  transition_probability_matrix_K80(t,alpha,beta):
#     q = np.zeros((4, 4))
#     p = np.zeros((4, 4))
#     q[:][:] = beta
#     q[0][1] = q[1][0] = q[2][3] = q[3][2] = alpha
#     q[0][0] = q[1][1] = q[2][2] = q[3][3] = -(alpha + 2 * beta)
#     # pprint.pprint(q)
#     eigval, eigvec = la.eig(q)
#     eigvec_inv = la.inv(eigvec)
#     p = np.dot(eigvec,np.dot(np.diag(np.exp(eigval * t)), eigvec_inv))
#     return p
#===================================================================================================================
# v = np.arange( 0 , 1.1, 0.01)
# y1 = []
# y2 = []
# y3 = []
# d = []
# for i in v:
#     p = transition_probability_matrix_K80(i,i,2*i)
#     y1.append(p[0][0])
#     y2.append(p[0][1])
#     y3.append(p[0][2])
#     d.append((p[0][1]+2*p[0][2])* i)
#
#
#
# pprint.pprint(p)
# plt.plot(d, y1, label='P(AA)')
# plt.plot(d, y2, label='P(AC)')
# plt.plot(d, y3, label='P(AG)')
# plt.ylabel('probablity')
# plt.xlabel('v')
# plt.axhline(q[0][0], color='r', ls='-.')
# plt.legend(loc='lower right')
# plt.show()


# alignment = dendropy.DnaCharacterMatrix.get\
#     (file=open("/home/nehleh/0_Research/PhD/Data/test.fasta"), schema="fasta" , )
# print(alignment.sequence_size)
#
# v = np.arange( 0 , 2, 0.01)
# y1 = []
# y2 = []
# d = []
# for i in v:
#     # p = transition_probability_matrix_JC(i,i)
#     p = transition_probability_matrix_GTR(alignment,i)
#     y1.append(p[0][0])
#     y2.append(p[0][1])
#     d.append(i * i)
#
#
#
# pprint.pprint(p)
# plt.plot(d, y1, label='P(AA)')
# plt.plot(d, y2, label='P(AC)')
# plt.ylabel('probablity')
# plt.xlabel('v')
# plt.axhline(p[0][0], color='r', ls='-.')
# plt.legend(loc='lower right')
# plt.show()



# alignment = ['AAAAAGGCAA','GGGCTCTTAA']
# r = ratio_matrix_computation(alignment)
# n = len(alignment[0])
# pprint.pprint(r)
# print(r['A']['G'])
# print(avg_base_frequencies(['ACCT','ACGT']))

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

