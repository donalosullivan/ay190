from numpy import array,matrix,pi,exp,sin,linspace,log10
from numpy.fft import fft
import matplotlib.pyplot as plt
from mycode import dft

i=1j
def f(x): return sin(x)

###PART 1
X = linspace(0.01,2*pi,100)
f_X = f(X)
dft_X = dft(f_X)
fft_X = fft(f_X)

fig = plt.figure() 
ax1 = fig.add_subplot(311)
ax1.plot(dft_X,label="My DFT")
ax1.legend()
ax2 = fig.add_subplot(312)
ax2.plot(fft_X,label="Numpy FFT")
ax2.legend()
ax3 = fig.add_subplot(313)
ax3.plot(dft_X-fft_X,label="Residual")
ax3.legend()
plt.savefig("fft_vs_dft.png")




