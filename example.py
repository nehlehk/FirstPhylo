import numpy


with open('/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/seq-gen-order/test', 'r') as f:
    data = f.read().replace('[', ' ')
    data = data.split(",")
    ll = []
    for elem in data:
        ll.append(float(elem))

signal = numpy.array(ll)

print(signal[7])


