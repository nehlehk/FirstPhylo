import os


help = "/home/nehleh/0_Research/Software/fastsimbac/fastSimBac_linux/example_input"

fastsimbac = "/home/nehleh/anaconda3/pkgs/fastsimbac-1.0.1_bd3ad13d8f79-h6dcb523_0/bin/fastSimBac "


param = "10 5000 -T -t .001 -r .00001 500 > fastsimbactree.txt" \
#         " -R /home/nehleh/0_Research/Software/fastsimbac/fastSimBac_linux/example_input/hotspot.txt" \
#         " -F /home/nehleh/0_Research/Software/fastsimbac/fastSimBac_linux/example_input/ascertainment.txt 0 2>trees.txt" \
#          " | /home/nehleh/anaconda3/pkgs/fastsimbac-1.0.1_bd3ad13d8f79-h6dcb523_0/bin/msformatter > 000haplotypes.txt"


cmd = fastsimbac + param
print(cmd)
os.system(cmd)