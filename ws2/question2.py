import numpy as np
from matplotlib import pyplot as plt

def f(x): return x**3 - 5*x**2 + x

def df(x): return 3*x**2 - 10*x +1

def df_for(x,h): return (f(x+h)-f(x))/h

def df_cen(x,h): return (f(x+h)-f(x-h))/(2*h)

h1 = np.float32(0.1)
h2 = h1/2

x1 = np.arange(-2,6,h1)
x2 = np.arange(-2,6,h2)
actual1 = df(x1)
actual2 = df(x2)

abserr1for = abs(df_for(x1,h1) - actual1)
abserr2for = abs(df_for(x2,h2) - actual2)

abserr1cen = abs(df_cen(x1,h1) - actual1)
abserr2cen = abs(df_cen(x2,h2) - actual2)

fs = 24

plt.figure()
plt.title("Forward Differencing Method")
plt.plot(x1,abserr1for,label="h=%.3f" % h1)
plt.plot(x2,abserr2for,label="h=%.3f" % h2)
plt.ylabel("Absolute error",fontsize=fs)
plt.xlabel("x",fontsize=fs)
plt.legend(loc=9)
plt.savefig("question2fig1.png")


plt.figure()
plt.title("Central Differencing Method")
plt.plot(x1,abserr1cen,label="h=%.3f" % h1)
plt.plot(x2,abserr2cen,label="h=%.3f" % h2)
plt.ylabel("Absolute error",fontsize=fs)
plt.xlabel("x",fontsize=fs)
plt.legend(loc=7)
plt.savefig("question2fig2.png")
