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

def upwind(dx,dt,v,y):
    y2 = y
    y2[1:] = y[1:] - v*(dt/dx)*(y[1:] - y[:-1])
    y2[0] = y2[1] #outflow boundary condition
    return y2     

def ftcs(dx,dt,v,y):
    y2 = y
    y2[1:-1] = y[1:-1] - ((v*dt)/(2*dx))*( y[2:] - y[:-2] )
    y2[0] = y2[1] #outflow boundary conditions
    y2[-1] = y2[-2]
    return y2
    
x = np.linspace(0,100,500) #Set up spatial range
x2 = np.linspace(0,100,75)
dx = x[1]-x[0] #spatial step size
v = 0.1 #advection velocity

n = len(x) #store range length
y = np.zeros(n) #initialize function as zeros
cfl = 1.0 #Courant-Friedrics-Levy factor
dt = cfl*dx/v #temporal step size
t = 0.0 #start t at 0

sigma = np.sqrt(15.0) #Gaussian parameters given in question
x0 = 30.0

print v*(dt/dx)
y = analytic(x0,x,v,t,sigma) #Set up initial conditions

# evolve (and show evolution)
mpl.ion()
mpl.figure()
mpl.plot(x,y,'x-') # numerical data
mpl.plot(x,analytic(x0,x,v,t,sigma),'r-') # analytic data
mpl.show()

yold2 = y #Store previous previous step
yold = y #Store previous step
ntmax = 1000
T,ERR = [],[]
for it in range(ntmax):

    t += dt #increment time

    yold2 = yold #shift data back one 'register'
    yold = y

    # get new data; ideally just call a function
    y = ftcs(dx,dt,v,y)

    # after update, apply boundary conditions
    # apply_bcs(x,y) 
    #[FILL IN CODE]

    # get analytic result for time t
    yana = analytic(x0,x,v,t,sigma)
    yana2 = analytic(x0,x2,v,t,sigma)
    # compute error estimate
    err = 0
    err = np.sum(abs(y-yana))
    T.append(t)
    ERR.append(np.log10(err))
    
    #COMMENT THIS OUT TO SWITCH TO OTHER PLOT MODE
    #mpl.clf()
    #mpl.plot(x,y,'b-')
    #mpl.plot(x2,yana2,'ro')
    #mpl.draw()
  
    #UN-COMMENT TO SWITCH PLOT MODE
    if it == 200:
        mpl.subplot(313)
        mpl.plot(x,y,'b-',label= "t = %.2f" % t )
        mpl.plot(x2,yana2,'ro')
        mpl.legend()
        mpl.xlabel("x",fontsize=20)
        mpl.draw()
        mpl.savefig("timeplot_ftcs.png")    
    elif it == 100:
        mpl.subplot(312)
        mpl.plot(x,y,'b-',label= "t = %.2f" % t )
        mpl.plot(x2,yana2,'ro')
        mpl.legend()
        mpl.draw()
    elif it==10:
        mpl.clf()
        mpl.subplot(311)
        mpl.plot(x,y,'b-',label= "t = %.2f" % t )
        mpl.plot(x2,yana2,'ro')
        mpl.legend()
        mpl.draw()
  

               
fs = 20
mpl.figure()
mpl.plot(np.log10(T),ERR)
mpl.xlabel(r"$Log_{10}(t)$",fontsize=fs)
mpl.ylabel(r"$Log_{10}( Err )$",fontsize=fs)
mpl.savefig("errplot_ftcs.png")
mpl.show()
