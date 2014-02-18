import sys,math
import numpy as np
import matplotlib.pyplot as mpl
import scipy as sp

def apply_bcs(x,y):
    # apply boundary conditions
    # you need to fill in code
    #[ FILL IN CODE ]
    return y


def analytic(x0,x,v,t,sig): return np.exp( -((x-x0-v*t)**2)/(2*sig**2))

    

# set up the grid here. Use a decent number of zones;
# perhaps to get a dx of 0.1
x = np.linspace(0,100,1000)
# parameters
dx = x[1]-x[0]
v = 0.1

n = len(x)
y = np.zeros(n)
cfl = 1.0
dt = cfl*dx/v
t = 0.0

# for initial data
sigma = np.sqrt(15.0)
x0 = 30.0

#set up initial conditions
y = analytic(x0,x,v,t,sigma)

# evolve (and show evolution)
mpl.ion()
mpl.figure()
mpl.plot(x,y,'x-') # numerical data
mpl.plot(x,analytic(x0,x,v,t,sigma),'r-') # analytic data
mpl.show()

yold2 = y
yold = y
ntmax = 2000
for it in range(ntmax):
    t += dt
    # save previous and previous previous data
    yold2 = yold
    yold = y

    # get new data; ideally just call a function
    #y = ????
    #[FILL IN CODE]

    # after update, apply boundary conditions
    # apply_bcs(x,y) 
    #[FILL IN CODE]

    # get analytic result for time t
    yana = analytic(x0,x,v,t,sigma)
    # compute error estimage
    err = 0
    # err = ???
    #[FILL IN CODE]
    print "it = ",it,err
    mpl.clf()
    # plot numerical result
    #[FILL IN CODE]
    #plot analytic results
    mpl.plot(x,yana,'r-')
    mpl.draw()


mpl.show()


