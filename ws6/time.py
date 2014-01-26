import matplotlib.pyplot as plt
from numpy import arange,log10,array
from timeit import timeit
from scipy.optimize import curve_fit as fit

def line(x,m,c): return m*x + c

###PART 2
N = arange(10,100,1)
tdft,tfft = [],[]
for n in N:
    tdft.append(timeit("dft(x)", number=10, \
                     setup="from mycode import dft;\
                     import pylab; x=pylab.rand(%d)" % n))
    tfft.append(timeit("fft(x)", number=10, \
                     setup="from numpy.fft import fft;\
                     import pylab; x=pylab.rand(%d)" % n))

#Logged data easier to compare power laws              
logN = array(log10(N))
logdft = log10(tdft)
logfft = log10(tfft)

#Fit line to DFT times to get slope
m,c = fit( line, logN, logdft )[0]
model = line(logN,m,c)

#Plot             
fs = 20
plt.figure()
plt.plot(logN,logdft,'bo',label="My DFT")
plt.plot(logN,logfft,'ro',label="FFT")
plt.plot(logN,model,'k-',label="y = %.2fx + %.2f" % (m,c))
plt.xlabel(r"$Log_{10}(N)$",fontsize=fs)
plt.ylabel(r"$Log_{10}(t)$",fontsize=fs)
plt.legend(loc=5)
plt.savefig("timing.png")
