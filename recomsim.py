import os
import numpy as np


# =====================================  fastsimbactree ====================================
seq_len = []
seq_tree = []
f = open("/home/nehleh/PycharmProjects/FirstPhylo/fastsimbactree.txt", "r")
treefile= "/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/5000/tree.tree"
for line in f:
    if  line.find('[') == 0:
        info = line
        tree = info.split(']')
        seq_len.append(tree[0].strip('['))
        seq_tree.append(tree[1])
        fi = open(treefile, "w")
        fi.write(tree[1])


# =====================================  seqgen ====================================
for i in np.arange(len(seq_len)):
    model = 'GTR'
    seqnum = seq_len[i]
    frequencies = '0.2184,0.2606,0.3265,0.1946'
    rates = '2.0431,0.0821,0,0.067,0,1'
    # frequencies = '0.25,0.25,0.25,0.25'
    # rates = '1,1,1,1,1,1'
    outfile = '/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/5000/'+str(i)+'.fasta'
    seqgen = '/home/nehleh/0_Research/Software/Seq-Gen-1.3.4/source/seq-gen'
    cmd = seqgen + '  -m'+model+ '  -l'+ seqnum + '  -f'+frequencies + '  -r'+rates +'  -of'+ '  <'+treefile+'>  ' + outfile
    print(cmd)
    os.system(cmd)

# # =====================================  combine fasta files ====================================
#
#
#     castfasta = '/home/nehleh/0_Research/Software/catfasta2phyml-master/catfasta2phyml.pl '
#     paramcast = ' -f  /home/nehleh/0_Research/PhD/Data/simulationdata/recombination/5000/*.fasta > /home/nehleh/0_Research/PhD/Data/simulationdata/recombination/5000/5000-recombin.fasta'
#     cmdcast = castfasta + paramcast
#     print(cmdcast)
#     os.system(cmdcast)
