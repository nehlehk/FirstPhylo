import dendropy
from dendropy.simulate import treesim


tips_number = 4

temp = []
for i in range(tips_number):
    temp.append(str(i+1))


taxa = dendropy.TaxonNamespace(temp)
tree = treesim.birth_death_tree(birth_rate=1.0 , death_rate= 0.5 ,ntax= tips_number , taxon_namespace = taxa)
tree.print_plot()
print(tree.as_string(schema='newick'))
f = open("/home/nehleh/0_Research/PhD/Data/tree/simtree1.tree", "w")
# "I planted a tree here :)"
f.write(tree.as_string(schema='newick'))
f.close()



