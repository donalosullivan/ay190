import numpy.random as random
import matplotlib.pyplot as plt
import numpy as np

def f(x): return x**2 + 1

def I_f(x): return (x**3)/3 + x
X = np.linspace(2,3,1000)
Y = f(X)

a,b = X[0],X[-1]
fa,fb = Y[0],Y[-1]

A1 = (b-a)*fb
A2 = (b-a)*(fb-fa)
Nvals = np.linspace(1e2,1e4,1000)
Is = []
for N in Nvals:
    n = 0
    for x in range(int(N)):

        xi = a + (b-a)*random.rand()
        yi = fa + (fb-fa)*random.rand()
        if yi <= f(xi): n+=1

    I_mc = A2*(float(n)/N) + (A1-A2)
    Is.append(I_mc)

I_an = I_f(b) - I_f(a)
Is = np.array(Is)
I_err = abs(Is - I_an)/I_an
fs = 20
plt.figure()
plt.plot(Nvals,I_err,'ko')
plt.xlabel("N",fontsize=fs)
plt.ylabel("% Error",fontsize=fs)
plt.savefig("q3fig.png")
plt.show()


