#!/usr/bin/env python

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from numpy import log10 as log
fs=20
def dudx(x): return 12*x - 4

def calc_rhs(u,x):
    # RHS Routine
    rhs = np.zeros(2)
    rhs[0] = u
    rhs[1] = dudx(x)
    return rhs

def integrate_FE(A,z,x,dx,res):
    # Forward-Euler integrator
    yy = np.zeros((res,2))
    yy[0,0] = A # boundary value A for y at x=0
    yy[0,1] = z # guessed boundary value for y' at x=0
    for i in range(res-1):
        yy[i+1,:] = yy[i,:] + dx*calc_rhs(yy[i,1],x[i])
    return yy
    
def rk_RHS(u,x):
    rhs = np.zeros(2)       
    rhs[0] = u
    rhs[1] = 12*x-4   
    return rhs
    
def rk2_integrate(A,z,x,dx,res):
    # Runge-Kutta
    yy = np.zeros((res,2)) 
    yy[0,0] = A # boundary value A for y at x=0
    yy[0,1] = z # guessed boundary value for y' at x=0
   
    for i in range(res-1):
        
        k1 = dx*rk_RHS(yy[i,1], x[i])
        k2 = dx*rk_RHS(yy[i,1]+k1[1]/2, x[i] + dx/2)
        yy[i+1,:] = yy[i,:] + k2
    return yy

def y(x): return 2.0*x**3 - 2*x**2 + 0.1*x
def bvp(integrator,lab="Unknown",col="r"):

    # set up grid
    xmin = 0.0
    xmax = 1.0
    npoints = 100
    x = np.linspace(xmin,xmax,npoints) # set up grid
    dx = x[1]-x[0] # step size

    # boundary values
    A = 0 # inner boundary
    B = 0.1 # outer boundary

    # get initial guess for derivative
    z0 = -1100000.0
    z1 = 10000000.0
    yy0 = integrator(A,z0,x,dx,npoints)
    yy1 = integrator(A,z1,x,dx,npoints)
    phi0 = yy0[npoints-1,0] - B
    phi1 = yy1[npoints-1,0] - B
    dphidz = (phi1-phi0)/(z1-z0)# dphi/dz

    i = 0
    itmax = 100
    err = 1.0e99 #Placeholder for error of each iteration
    criterion = 1.0e-15 #Desired accuracy
    errs = []
    z0 = z1
    phi0 = phi1
    while (err > criterion and i < itmax):
        z1 = z0 - phi0/dphidz # secand update
        yy = integrator(A,z1,x,dx,npoints)
        phi1 = yy[npoints-1,0] - B
        dphidz = (phi1-phi0)/(z1-z0) # dphi/dz numerical
        err = np.abs((z1-z0)/z0) # your error measure
        errs.append(err)
        z0 = z1
        phi0 = phi1
        i = i+1
    its = np.arange(1,i+1)
    plt.subplot(212)
    plt.plot(its,log(errs),"%so" % col,label=lab)
    plt.plot(its,log(errs),"%s-" % col)
    plt.xlabel(r"$Iteration$",fontsize=fs)
    plt.ylabel(r"$Log_{10}(Error)$",fontsize=fs)
    plt.legend()
    plt.subplot(211)
    plt.plot(x,yy[:,0],col,label=lab)


plt.figure()
bvp(integrate_FE,lab="FE",col="r")
bvp(rk2_integrate,lab="RK2",col="b")
xmin = 0.0
xmax = 1.0
x = np.linspace(xmin,xmax,25)
plt.plot(x,y(x),"ko",label = "Analytical")
plt.xlabel(r"$x$",fontsize=fs)
plt.ylabel(r"$y$",fontsize=fs)
plt.legend()
plt.show()




