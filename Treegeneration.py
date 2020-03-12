
import dendropy
from dendropy.simulate import treesim


tips_number = 20

temp = []
for i in range(tips_number):
    temp.append(str(i+1))


taxa = dendropy.TaxonNamespace(temp)
t = treesim.birth_death_tree(birth_rate=1.0 , death_rate= 0.5 ,ntax= tips_number , taxon_namespace = taxa)
t.print_plot()
print(t)





