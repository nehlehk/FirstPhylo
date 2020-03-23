# import libraries
import numpy
import ruptures as rpt
import matplotlib.pyplot as plt





with open('/home/nehleh/0_Research/PhD/Data/simulationdata/likelihoos_JC', 'r') as f:
    data = f.read().split(" ")
    ll = []
    for elem in data:
        try:
            ll.append(float(elem))
        except ValueError:
            pass
    print(ll[0:10])


signal = numpy.array(ll)

print(signal.shape)


# change point detection
model = "l2"   # "l1", "rbf", "linear", "normal", "ar"
# my_bkps = rpt.Dynp(model=model, min_size=3).fit_predict(signal,n_bkps=3)
my_bkps = rpt.Window(model=model, width= 15).fit_predict(signal,n_bkps=5)
print(my_bkps)
# show results
rpt.show.display(signal, my_bkps , figsize =(8,3))
plt.show()
