import os


bootstrapp = " -m GTRCAT -p 12345 -b 12345  -V -# 100 -s /home/nehleh/0_Research/PhD/Data/LL_vector/JC69_100.fasta  -n T14"
bipartition = " -m GTRCAT -p 12345 -f b -t RAxML_bestTree.T13 -z RAxML_bootstrap.T14 -n T15"

help = "/home/nehleh/0_Research/Software/standard-RAxML-master/raxmlHPC -h"

raxml = "/home/nehleh/0_Research/Software/standard-RAxML-master/raxmlHPC"

param = " -m GTRGAMMA   -p 12345 -s /home/nehleh/0_Research/PhD/Data/Ecoli-10_Aligned/10_aligned.fasta -N 10 -n RaxMLtree"

param2 = " -m GTRGAMMA  --JC69  -p 12345 -s /home/nehleh/0_Research/PhD/Data/Ecoli-10_Aligned/10_aligned_reduced.fasta -n likelihood_JC "

paramlikelihood = " -m GTRGAMMA   -p 12345 -s /home/nehleh/0_Research/PhD/Data/Ecoli-10_Aligned/10_aligned_reduced.fasta -n likelihood_GTR"

likelihood = " -f g -z /home/nehleh/0_Research/PhD/Data/Ecoli-10_Aligned/RAxML_bestTree.RaxMLtree"

cmd = raxml + param2 + likelihood
print(cmd)
os.system(cmd)