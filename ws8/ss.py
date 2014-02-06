#!/usr/bin/env python

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from numpy import log10 as log
# global constants
ggrav = 6.67e-8
msun  = 1.99e33

# EOS parameters
# for white dwarfs:
polyG = 4.0/3.0
polyK = 1.244e15*0.5**polyG


##DIFFERENTIAL EQS
def y_press(m,rho,r): return -ggrav*m*rho/r**2
def y_mass(rho,r): return 4*np.pi*rho*r**2

#######################################
# function definitions
def tov_RHS(rad,p,rho,m):
    
    # RHS function   
    rhs = np.zeros(2)
    if(rad > 1.0e-10):
        rhs[0] = y_press(m,rho,rad) #FILL3
        rhs[1] = y_mass(rho,rad) #FILL4
    else:
        rhs[0] = 0.0
        rhs[1] = 0.0

    return rhs

def tov_integrate_FE(rad,dr,p,rho,m):

    # Forward-Euler Integrator

    new = np.zeros(2) # New vals for [P,M]
    old = np.zeros(2) # Old vals for [P,M]
    old[0] = p
    old[1] = m

    # forward Euler integrator
    new = old + dr*tov_RHS(rad,p,rho,m) #FILL5
    
    # assign outputs
    pnew = new[0]
    mnew = new[1]
    
    return (pnew,mnew)

def rk_RHS(r,dr,p,rho,m):

    rhs = np.zeros(2)
    
    if r < 1.0e-10: return rhs
    
    #K1 values
    k1_p = dr*y_press(m,rho,r) #Pressure
    k1_m = dr*y_mass(rho,r) #Mass
    #print "%10.3e %10.3e %10.3e %10.3e" %(k1_p, m,rho,r)
    
    p2 = p + k1_p/2 #Get intermediate value of p
    rho2 = (abs(p2/polyK))**(1/polyG) #Get corresponding rho
    
    #K2 values
    k2_p = dr*y_press(m + k1_m/2, rho2, r + dr/2)
    k2_m = dr*y_mass(rho2, r + dr/2)
    
    rhs[0] = k2_p
    rhs[1] = k2_m
    
    return rhs
    
def rk2_integrate(rad,dr,p,rho,m):
       
    new,old = np.zeros(2),np.zeros(2)
    old[0],old[1] = p,m
      
    new = old + rk_RHS(rad,dr,p,rho,m)
    pnew = new[0]
    mnew = new[1]
    
    return (pnew,mnew)
    
#######################################

# set up grid
npoints = 1000

#Returns [r,rho(r),P(r),M(r),nsurf] for the white dwarf
def wd(integrator,npoints):

    radmax = 2.0e8 # 2000 km
    radius = np.linspace(0,radmax,npoints) #FILL1
    dr = radius[1]-radius[0]

    # set up variables
    press = np.zeros(npoints)
    rho   = np.zeros(npoints)
    mass  = np.zeros(npoints)

    # set up central values
    rho[0]   = 1.0e10
    press[0] = polyK * rho[0]**polyG
    mass[0]  = 0.0

    # set up termination criterion
    press_min = 1.0e-10 * press[0]

    nsurf = 0
    for n in range(npoints-1):
        
        (press[n+1],mass[n+1]) = integrator(radius[n],
                                                  dr,
                                                  press[n],
                                                  rho[n],mass[n])
        # check for termination criterion
        if(press[n+1] < press_min and nsurf==0):
            nsurf = n

        if(n+1 > nsurf and nsurf > 0):
            press[n+1] = press[nsurf]
            rho[n+1]   = rho[nsurf]
            mass[n+1]  = mass[nsurf]

        # invert the EOS to get density
        rho[n+1] = (press[n+1]/polyK)**(1/polyG) #FILL2


    return [radius,rho,press,mass,nsurf]

fs = 24

def part2():
    R,rho,P,M,nsurf = wd(tov_integrate_FE,1000)
    print "Mass: %10.5e\nRadius: %10.5e" % (M[nsurf]/msun,R[nsurf]/1e5)
    
#PART 3 - PLOTTING CONVERGENCE
def part3():
    eul_vals = []
    rk2_vals = []
    npoints_range = np.arange(200,10000,100)
    for res in npoints_range:
        print ("%10i" % res),
        
        R,rho,P,M,nsurf = wd(tov_integrate_FE,res)
        m = M[nsurf]/msun
        eul_vals.append(m)
        print "%10s %10.3e" % ("EUL",m),
        
        R,rho,P,M,nsurf = wd(rk2_integrate,res)
        m = M[nsurf]/msun
        rk2_vals.append(m)
        print "%10s %10.3e" % ("RK2",m)
    eul_vals = np.array(eul_vals)
    rk2_vals = np.array(rk2_vals)
    #npoints_range = log(npoints_range)


    plt.figure()
    plt.plot(npoints_range,rk2_vals,'b',label="RK2")
    plt.plot(npoints_range,eul_vals,'r',label="Euler")
    plt.xlabel(r"$Resolution$",fontsize=fs)
    plt.ylabel(r"$Mass (M_{\odot})$",fontsize=fs)
    plt.legend(loc=5)
    plt.savefig("convergence.png")
    plt.show()

#PART 4 - PLOTTING PROFILES
def part4():
    R,rho,P,M,nsurf = wd(rk2_integrate,10000)

    Mnorm = M[1:nsurf]/M[nsurf]
    Pnorm = log(P[1:nsurf]/P[0])
    rhonorm = log(rho[1:nsurf]/rho[0])
    R = R[1:nsurf]/1e5
    plt.figure()
    plt.plot(R,Pnorm,'b',label=r"$f = P / P_{max}$")
    plt.plot(R,rhonorm,'r',label=r"$f = \rho / \rho_{max}$")
    plt.xlabel(r"$R[km]$",fontsize=fs)
    plt.ylabel(r"$Log_{10}(f)$",fontsize=fs)
    plt.legend(loc=8)
    y2 = plt.twinx()
    plt.ylabel(r"$g$",fontsize=fs)
    plt.plot(R,Mnorm,'k',label=r"$g = M / M_{tot}$")
    plt.legend()
    plt.savefig("profiles.png")
    plt.show()


part2()



