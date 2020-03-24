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


signal = numpy.array(ll)




# change point detection
model = "l1"   # "l1", "l2", "rbf", "linear", "normal", "ar"

# my_bkps = rpt.Dynp(model=model, min_size=3).fit_predict(signal,n_bkps=3)

# my_bkps = rpt.Window(model=model, width= 15).fit_predict(signal,n_bkps=5)

my_bkps = rpt.Pelt(model = model, min_size=5, jump=100).fit_predict(signal,pen=3*numpy.log(signal.shape[0]))



print(my_bkps)
# show results
rpt.show.display(signal, my_bkps , figsize =(15,7))
plt.show()
