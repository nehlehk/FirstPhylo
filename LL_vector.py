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

tree = Tree.get_from_path('/home/nehleh/0_Research/PhD/Data/LL_vector/tree.tree', 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/0_Research/PhD/Data/LL_vector/JC69_100.fasta"), schema="fasta")


tips = len(alignment)
alignment_len = alignment.sequence_size

dna = []
for i in range(tips):
    dna.append(alignment[i][0])



# partials = numpy.zeros((2 * tips - 1, alignment_len, 4))
# print(partials.shape)



partial = numpy.zeros((4,2*tips))
# print(partial.shape)


for i in range(4):
    for j in range(4):
        p[i][j] = 0.25 - 0.25* math.exp(-4*br_length/3)
    p[i][i] = 0.25 + 0.75*math.exp(-4*br_length/3)
# pprint.pprint(p)


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
        print(dna[pos])
        pos += 1

        # i = give_index(dna[node.index-1])
        # partial[i][node.index-1] = 1
    # else:
    #         for j in range(4):
    #             sump = []
    #             for x in node.child_node_iter():
    #                 z = 0
    #                 for k in range(4):
    #                    z  += p[j][k] * partial[k][x.index-1]
    #                 sump.append(z)
    #             partial[j][node.index-1] = numpy.prod(sump)


# print(partial)
#
# ll = sum(partial[:,2*tips -3] * 0.25)
#
# print("likelihood = {} and log-likelihood = {} ".format(round(ll,7) , round(numpy.log(ll),7)))