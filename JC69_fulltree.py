import random
import dendropy




mle = dendropy.Tree.get_from_path('/home/nehleh/0_Research/PhD/Data/pythonidae.mle.nex', 'nexus')

mytree = dendropy.Tree.get_from_path('/home/nehleh/0_Research/PhD/Data/tree.tre', 'newick')

mytree2 = dendropy.Tree.get_from_path('/home/nehleh/0_Research/PhD/Data/treelength.tre', 'newick')


# print(mle.length())

mrca = mytree.mrca(taxon_labels = ['node1' ,'node2'])

print(mrca.taxon)


for node in mytree.postorder_node_iter():
    print(node.taxon)


print(len(mytree.seed_node.adjacent_nodes()))



leaves= mytree.leaf_nodes()

for  i, lf in enumerate(leaves):
    print(leaves[i].taxon)



mytree.print_plot()

#  checking for sth







