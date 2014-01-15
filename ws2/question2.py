import numpy as np
from matplotlib import pyplot as plt

def f(x): return x**3 - 5*x**2 + x

def df(x): return 3*x**2 - 10*x

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


plt.figure()
plt.plot(x1,abserr1for,label="for h1")
plt.plot(x2,abserr2for,label="for h2")
plt.legend()
plt.show()


plt.figure()
plt.plot(x1,abserr1cen,label="cen h1")
plt.plot(x2,abserr2cen,label="cen h2")
plt.legend()
plt.show()	

