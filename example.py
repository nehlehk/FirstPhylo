import os
treefile= "/home/nehleh/0_Research/PhD/Data/simulationdata/recombinant.tree"

seqnum = '2000'
model = 'GTR'
partion = '5'
frequencies = '0.2184,0.2606,0.3265,0.1946'
rates = '2.0431,0.0821,0,0.067,0,1'
# frequencies = '0.25,0.25,0.25,0.25'
# rates = '1,1,1,1,1,1'
outfile = '/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/5000/recom_test.fasta'
seqgen = '/home/nehleh/0_Research/Software/Seq-Gen-1.3.4/source/seq-gen'
cmd = seqgen + '  -m'+model+ '  -l'+ seqnum + '  -f'+frequencies  + '  -p' +  partion  + '  -r'+rates +'  -of'+ '  <'+treefile+'>  ' + outfile
print(cmd)
os.system(cmd)