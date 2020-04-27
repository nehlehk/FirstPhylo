import os
import numpy as np

seq_len = []
seq_tree = []
f = open("/home/nehleh/PycharmProjects/FirstPhylo/fastsimbactree.txt", "r")
for line in f:
    if  line.find('[') == 0:
        info = line
        tree = info.split(']')
        seq_len.append(tree[0].strip('['))
        seq_tree.append(tree[1])




for i in np.arange(len(seq_len)):
    model = 'GTR'
    seqnum = seq_len[0]
    frequencies = '0.2184,0.2606,0.3265,0.1946'
    rates = '2.0431,0.0821,0,0.067,0,1'
    # frequencies = '0.25,0.25,0.25,0.25'
    # rates = '1,1,1,1,1,1'
    outfile = '/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/5000/'+str(0)+'.fasta'
    seqgen = '/home/nehleh/0_Research/Software/Seq-Gen-1.3.4/source/seq-gen'
    cmd = seqgen + '  -m'+model+ '  -l'+ seqnum + '  -f'+frequencies + '  -r'+rates +'  -of'+ '  <'+seq_tree[0]+'>  ' + outfile
    print(cmd)
    os.system(cmd)



/home/nehleh/0_Research/Software/Seq-Gen-1.3.4/source/seq-gen  -mGTR  -l1287  -f0.2184,0.2606,0.3265,0.1946  -r2.0431,0.0821,0,0.067,0,1  -of  <(((6:0.00442254,0:0.00442254):0.180659,(3:0.0219039,(2:0.00972834,1:0.00972834):0.0121755):0.163178):1.0284,(8:0.14008,(9:0.0883774,((5:0.0332952,4:0.0332952):0.000313537,7:0.0336087):0.0547686):0.0517026):1.0734);
>  /home/nehleh/0_Research/PhD/Data/simulationdata/recombination/5000/0.fasta


