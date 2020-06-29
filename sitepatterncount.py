from dendropy import Tree, DnaCharacterMatrix
import numpy
import math
import dendropy


alignment = dendropy.DnaCharacterMatrix.get(file=open("/media/nehleh/295eaca0-f110-4b79-8cbe-bc99f9f61cbd/nehleh/0_Research/PhD/Data/simulationdata/recombination/clonalframe/wholegenome.fasta"), schema="fasta")
alignment_len = alignment.sequence_size
tips = len(alignment)


coulmns = []
for l in range(alignment_len):
    col = ""
    for t in range(tips):
        col += str(alignment[t][l])
    coulmns.append(col)

myset = list(set(coulmns))
# print((myset))
print(len(myset))
# print(myset[0:10])
print("----------------------")
a = 0
for i in range(len(myset)):
    print(coulmns.count(myset[i]))
    a = a + coulmns.count(myset[i])

print(a)

# print(coulmns[0:10])
# print(coulmns.count(coulmns[6]))