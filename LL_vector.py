from dendropy import Tree, DnaCharacterMatrix
import numpy
import math
import dendropy
import pprint


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

br_length = 0.1
p = numpy.zeros((4, 4))
LL_vector = []

tree = Tree.get_from_path('/home/nehleh/0_Research/PhD/Data/LL_vector/test/tree.tree', 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/0_Research/PhD/Data/LL_vector/test/test.fasta"), schema="fasta")


tips = len(alignment)
alignment_len = alignment.sequence_size


for l in range(alignment_len):
    for i in range(tips):
        alignment[i][l]



# dna = 'GCAA'




partial = numpy.zeros((4, 2 * tips , alignment_len))


for i in range(4):
    for j in range(4):
        p[i][j] = 0.25 - 0.25* math.exp(-4*br_length/3)
    p[i][i] = 0.25 + 0.75*math.exp(-4*br_length/3)
# pprint.pprint(p)

for j in range(tips):
    col = ''
    for i in range(alignment_len):
        col += str(alignment[i][j])
        
        for node in tree.postorder_node_iter():
            node.index = -1
            node.annotations.add_bound_attribute("index")

        s = tips + 1
        for id,node in enumerate(tree.postorder_node_iter()):
            if not node.is_leaf():
                node.index = s
                s += 1
            else:
                for idx, name in enumerate(col):
                    if idx + 1 == int(node.taxon.label):
                        node.index = idx+1
                        break


    # pos = 0
    # for node in tree.postorder_node_iter():
    #     if node.is_leaf():
    #         i = give_index(str(dna[pos]))
    #         pos += 1
    #         # i = give_index(dna[node.index-1])
    #         partial[d][node.index][i] = 1
    #     else:
    #             for j in range(4):
    #                 sump = []
    #                 for x in node.child_node_iter():
    #                     z = 0
    #                     for k in range(4):
    #                        z  += p[j][k] * partial[d][x.index][k]
    #                     sump.append(z)
    #                 partial[d][node.index][j] = numpy.prod(sump)





#     ll = sum(partial[:,2*tips -1][d] * 0.25)
#
#     LL_vector.append(round(numpy.log(ll),7))
#
#     print("likelihood = {} and log-likelihood = {} ".format(round(ll,7) , round(numpy.log(ll),7)))
#

print(partial)
# print(LL_vector)