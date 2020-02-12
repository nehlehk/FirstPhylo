import dendropy
import numpy as np
import math
import numpy.linalg as la
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
column = 'CGCC'

s = 4 #tips number
br_length = 0.1
p = np.zeros((4, 4))

mytree = dendropy.Tree.get_from_path('/home/nehleh/0_Research/PhD/Data/tree.tre', 'newick')
partial = np.zeros((4,2*s -2))

for i in range(4):
    for j in range(4):
        p[i][j] = 0.25 - 0.25* math.exp(-4*br_length/3)
    p[i][i] = 0.25 + 0.75*math.exp(-4*br_length/3)

pprint.pprint(p)

post_order =[]


for  index, lf in enumerate(mytree.leaf_nodes()):
    i = give_index(column[index])
    partial[i][index] = 1



for index,node in enumerate(mytree.postorder_node_iter()):
    post_order.append(node.taxon)
    # if node.is_leaf():
    #     i = give_index(column[index])
    #     partial[i][index] = 1
    # else:
    if not node.is_leaf():
        children = node.child_nodes()
        leftchild =  post_order.index(children[0].taxon)
        rightchild = post_order.index(children[1].taxon)
        for j in range(4):
            left = 0
            right = 0
            for k in range(4):
                left += p[j][k] * partial[k][leftchild]
                right += p[j][k] * partial[k][rightchild]
            partial[j][index] = left * right



print(partial)





# mytree.print_plot()








