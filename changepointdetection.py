# import libraries
import numpy
import ruptures as rpt
import matplotlib.pyplot as plt




# creation of data
with open('/home/nehleh/0_Research/PhD/Data/simulationdata/likelihoos_JC', 'r') as f:
    data = f.read().split(" ")
    ll = []
    for elem in data:
        try:
            ll.append(float(elem))
        except ValueError:
            pass


signal = numpy.array(ll)

alignlen = 200000
mean = numpy.mean(signal)
std = numpy.std(signal)



# change point detection
model = "l1"   # "l1", "rbf", "linear", "normal", "ar"

# search_method = 'dynamic programming'
# my_bkps = rpt.Dynp(model=model, min_size=100).fit_predict(signal,n_bkps=5)


# search_method = 'Window-based change point detection'
# my_bkps = rpt.Window(model=model, width= 5000).fit_predict(signal,pen=numpy.log(alignlen)*mean*std**2)
# my_bkps = rpt.Window(model=model, width= 5000).fit_predict(signal,epsilon=alignlen*mean*std**2)
# my_bkps = rpt.Window(model=model, width= 10000).fit_predict(signal,epsilon=10)



search_method = 'Exact segmentation: Pelt'
my_bkps = rpt.Pelt(model = model, min_size=500).fit_predict(signal,pen=100)

print(my_bkps)



# show results

rpt.show.display(signal, my_bkps , figsize =(15,6))
plt.title(search_method , y = 0.1)
plt.show()


