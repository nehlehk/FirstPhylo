# import libraries
from Bio import AlignIO
from Bio.AlignIO import MauveIO
import pprint




align = AlignIO.parse("/home/nehleh/0_Research/PhD/Data/100k/GTR_100k", "mauve")

alignments = list(align)


for id in range(len(alignments)):
    # print("MYID************************",id)
    for idx,record in enumerate(alignments[id]):
        if record.name == '1':
            print(record)


lcb=arr = [[0 for i in range(5)] for j in range(len(alignments))]



for id in range(len(alignments)):
    # print("MYID************************",id)
    for idx,record in enumerate(alignments[id]):
            lcb[id][0] = record.name
            lcb[id][1] = record.annotations['start']
            lcb[id][2] = record.annotations['end']
            lcb[id][3] = record.seq
            lcb[id][4] = record.annotations['strand']


pprint.pprint(lcb)

name = []
for i in range(len(lcb)):
    name.append(lcb[i][0])

# print(name)
uniqename= list(set(name))

# for i in uniqename:
#     for id,item in enumerate(lcb):
#         if item[0] == i:
#             print(item)

#
# for record in alignments:
#         print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
#         print(record)