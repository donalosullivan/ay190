import sys,math
import numpy as np
import matplotlib.pyplot as mpl
import scipy as sp

def analytic(x,t,L): return np.sin( 2*np.pi*x/L )/8.0

def upwind(dx,dt,y):
    y2 = np.zeros(len(y))
    for i in range(1,len(y2)-1):
        if y[i]>=0:
            y2[i] = y[i] - y[i]*(dt/dx)*(y[i] - y[i-1])
        else:
            y2[i] = y[i] - y[i]*(dt/dx)*(y[i+1] - y[i])
    y2[0] = y2[1] #outflow boundary condition
    y2[-1] = y2[-2]
    return y2     
L = 100
x = np.linspace(0,L,500) #Set up spatial range
dx = x[1]-x[0] #spatial step size
n = len(x) #store range length
y = np.zeros(n) #initialize function as zeros

t = 0.0 #start t at 0
y = analytic(x,t,L) #Set up initial conditions
cfl = 1.0 #Courant-Friedrics-Levy factor
dt = cfl*dx/max(y) #temporal step size


# evolve (and show evolution)
mpl.ion()
mpl.figure()
mpl.plot(x,y,'x-') # numerical data
mpl.plot(x,analytic(x,t,L),'r-') # analytic data
mpl.show()

yold2 = y #Store previous previous step
yold = y #Store previous step
ntmax = 150
T,ERR = [],[]

for it in range(ntmax):

    t += dt #increment time

    yold2 = yold #shift data back one 'register'
    yold = y

    # get new data; ideally just call a function
    y = upwind(dx,dt,y)
    print it,dt
    # get analytic result for time t
    yana = analytic(x,t,L)

    # compute error estimate
    err = 0
    err = np.sum(abs(y-yana))
    T.append(t)
    ERR.append(np.log10(err))
    
    #COMMENT THIS OUT TO SWITCH TO OTHER PLOT MODE
    #mpl.clf()
    #mpl.plot(x,y,'b-')
    #mpl.draw()
  
    #UN-COMMENT TO SWITCH PLOT MODE
    if it == 90:
        mpl.subplot(313)
        mpl.plot(x,y,'b-',label= "t = %.2f" % t )
        mpl.legend()
        mpl.xlabel("x",fontsize=20)
        mpl.draw()
        mpl.savefig("timeplot_burger.png")    
    elif it == 60:
        mpl.subplot(312)
        mpl.plot(x,y,'b-',label= "t = %.2f" % t )
        mpl.legend()
        mpl.draw()
    elif it==10:
        mpl.clf()
        mpl.subplot(311)
        mpl.plot(x,y,'b-',label= "t = %.2f" % t )
        mpl.legend()
        mpl.draw()
  

