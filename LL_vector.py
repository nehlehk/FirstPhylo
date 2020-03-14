from dendropy import Tree, DnaCharacterMatrix
import numpy
import math
import dendropy
import numpy.linalg as la
import pprint

p = numpy.zeros((4, 4))


def give_index(c):
    if c == "A":
        return 0
    elif c == "C":
        return 1
    elif c == "G":
        return 2
    elif c == "T":
        return 3


# =======================================================================================================
def computelikelihood(tree,dna,model):

    partial = numpy.zeros((4, 2 * tips))

    if model == 'JC69':
        p = numpy.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                p[i][j] = 0.25 - 0.25 * math.exp(-4 * br_length / 3)
            p[i][i] = 0.25 + 0.75 * math.exp(-4 * br_length / 3)


    if model == 'GTR':
         mu = 0
         freq = numpy.zeros((4, 4))
         q = numpy.zeros((4, 4))
         sqrtPi = numpy.zeros((4, 4))
         sqrtPiInv = numpy.zeros((4, 4))
         exchang = numpy.zeros((4, 4))
         s = numpy.zeros((4, 4))
         fun = numpy.zeros(alignment_len)

         freq = numpy.diag(pi)
         sqrtPi = numpy.diag(numpy.sqrt(pi))
         sqrtPiInv = numpy.diag(1.0 / numpy.sqrt(pi))
         mu = 1 / (2 * ((a * pi[0] * pi[1]) + (b * pi[0] * pi[2]) + (c * pi[0] * pi[3]) + (d * pi[1] * pi[2]) + (
                 e * pi[1] * pi[3]) + (pi[2] * pi[3])))
         exchang[0][1] = exchang[1][0] = a
         exchang[0][2] = exchang[2][0] = b
         exchang[0][3] = exchang[3][0] = c
         exchang[1][2] = exchang[2][1] = d
         exchang[1][3] = exchang[3][1] = e
         exchang[2][3] = exchang[3][2] = f

         q = numpy.multiply(numpy.dot(exchang, freq), mu)

         for i in range(4):
             q[i][i] = -sum(q[i][0:4])

         s = numpy.dot(sqrtPi, numpy.dot(q, sqrtPiInv))

         eigval, eigvec = la.eig(s)
         eigvec_inv = la.inv(eigvec)

         left = numpy.dot(sqrtPi, eigvec)
         right = numpy.dot(eigvec_inv, sqrtPiInv)

         p = numpy.dot(left, numpy.dot(numpy.diag(numpy.exp(eigval * br_length)), right))



    for node in tree.postorder_node_iter():
        node.index = -1
        node.annotations.add_bound_attribute("index")

    s = tips + 1
    for id,node in enumerate(tree.postorder_node_iter()):
        if not node.is_leaf():
            node.index = s
            s += 1
        else:
            for idx, name in enumerate(dna):
                if idx + 1 == int(node.taxon.label):
                    node.index = idx+1
                    break
    pos = 0
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            i = give_index(str(dna[pos]))
            pos += 1
            # i = give_index(dna[node.index-1])
            partial[i][node.index] = 1
        else:
            for j in range(4):
                sump = []
                for x in node.child_node_iter():
                    z = 0
                    for k in range(4):
                       z  += p[j][k] * partial[k][x.index]
                    sump.append(z)
                partial[j][node.index] = numpy.prod(sump)

    # print(partial)

    ll = sum(partial[:,2*tips -1] * 0.25)

    # print("likelihood = {} and log-likelihood = {} ".format(round(ll,7) , round(numpy.log(ll),7)))
    return round(numpy.log(ll),7)
#=======================================================================================================================
rates = numpy.ones(5)
a, b, c, d, e = rates
f = 1
pi = [0.25]*4
br_length = 0.1
LL_vector_JC69 = []
LL_vector_GTR = []


tree = Tree.get_from_path('/home/nehleh/0_Research/PhD/Data/LL_vector/tree.tree', 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/0_Research/PhD/Data/LL_vector/JC69_100.fasta"), schema="fasta")


tips = len(alignment)
alignment_len = alignment.sequence_size

for l in range(alignment_len):
    col = ""
    for t in range(tips):
        col += str(alignment[t][l])
    LL_vector_JC69.append(computelikelihood(tree,col,'JC69'))
    LL_vector_GTR.append(computelikelihood(tree, col, 'GTR'))

print(LL_vector_JC69)
print(LL_vector_GTR)

