# import libraries
import numpy
import matplotlib.pyplot as plt
import ruptures as rpt



# creation of data
with open('/home/nkargarf/Desktop/likelihoos_JC', 'r') as f:
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
# my_bkps = rpt.Window(model=model, width= 15).fit_predict(signal,pen=numpy.log(alignlen)*mean*std**2)


# search_method = 'Exact segmentation: Pelt'
# my_bkps = rpt.Pelt(model = model, min_size=5, jump=100).fit_predict(signal,pen=3*numpy.log(signal.shape[0]))


# search_method = 'Bottom-up segmentation'
# my_bkps = rpt.BottomUp(model = model).fit_predict(signal,pen=numpy.log(alignlen)*mean*std**2)


search_method = 'Binary segmentation'
my_bkps = rpt.Binseg(model = model).fit_predict(signal,pen=100)



print(my_bkps)



# show results
rpt.show.display(signal, my_bkps , figsize =(15,7))
plt.title(search_method)
plt.show()
