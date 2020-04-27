import os

seq_count = '10'
seqnum = '5000'
# =====================================  fastsimbactree ====================================
help = "/home/nehleh/0_Research/Software/fastsimbac/fastSimBac_linux/example_input"

fastsimbac = "/home/nehleh/anaconda3/pkgs/fastsimbac-1.0.1_bd3ad13d8f79-h6dcb523_0/bin/fastSimBac  "

paramfast = seq_count + "  " + seqnum + " -T -t .001 -r .00001 500 > /home/nehleh/0_Research/PhD/Data/simulationdata/recombination/5000/fastsimbactree.txt"

cmdfastsimbac = fastsimbac + paramfast
print(cmdfastsimbac)
os.system(cmdfastsimbac)

# =====================================  prepare tree for seq-gen ====================================

f = open("/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/5000/fastsimbactree.txt", "r")
treefile= "/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/5000/tree.tree"
partition = 0
for line in f:
    if  line.find('[') == 0:
        partition += 1
        fi = open(treefile, "a")
        fi.write(line)
        fi.close()

# =====================================  seqgen ====================================
model = 'GTR'
frequencies = '0.2184,0.2606,0.3265,0.1946'
rates = '2.0431,0.0821,0,0.067,0,1'
# frequencies = '0.25,0.25,0.25,0.25'
# rates = '1,1,1,1,1,1'
outfile = '/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/5000/recombined.fasta'
seqgen = '/home/nehleh/0_Research/Software/Seq-Gen-1.3.4/source/seq-gen'
cmd = seqgen + '  -m'+model+ '  -l'+ seqnum + '  -f'+frequencies  + '  -p' + str(partition)  + '  -r'+rates +'  -of'+ '  <'+treefile+'>  ' + outfile
print(cmd)
os.system(cmd)



