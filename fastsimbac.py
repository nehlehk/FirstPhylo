import os


help = "/home/nehleh/0_Research/Software/fastsimbac/fastSimBac_linux/example_input"

fastsimbac = "/home/nehleh/anaconda3/pkgs/fastsimbac-1.0.1_bd3ad13d8f79-h6dcb523_0/bin/fastSimBac "


paramfast = "10 5000 -T -t .001 -r .00001 500 > /home/nehleh/PycharmProjects/FirstPhylo/fastsimbactree.txt"

cmdfastsimbac = fastsimbac + paramfast
print(cmdfastsimbac)
os.system(cmdfastsimbac)