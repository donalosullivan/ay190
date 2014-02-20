import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

ggrav = 6.67e-8
def phi(r,density,R): return (2.0/3.0)*np.pi*ggrav*density*(r**2 - 3*R**2)
#Part 1: Loading data

#1.1 Load and store
print "Loading Data..." 
data = np.loadtxt("presupernova.dat") #Parse data file
print "Done."
r = data[:,2] #Identified column 3 as r values because it's increasing, positive and starts <10^9
dens = data[:,3] #Identified column 4 as density because it's decreasing and positive
mass = data[:,1] #Mass, idenfitied by being correct order of magnitude

#1.2 Find index of outermost radius where r < 1e9
edge = 0
for i in range(len(r)):
    if r[i]>1e9: 
        edge = i-1
        break;
       
#1.3 Truncate data to r < 1e9
r = r[:edge]
dens = dens[:edge]
mass = mass[:edge]


#Part 2: Interpolating

#2.1 Interpolate
print "Interpolating data..."
dens_func = interp1d(r,dens,kind="cubic") #Scipy's interpolate using Piecewise Cubic Hermite
print "Done."

#2.2 Generate model
res = 1000
r_grid = np.linspace(r[0],r[-1],res) #Create evenly spaced grid
dens_grid = dens_func(r_grid) #interpolated density function

#2.3 Plot
fs = 20
plt.figure()
plt.plot(np.log10(r_grid),np.log10(dens_grid),'r.',label="Interpolated")
plt.plot(np.log10(r),np.log10(dens),'k-',label="Data")
plt.xlabel(r"$Log_{10}( R [cm] )$",fontsize=fs)
plt.ylabel(r"$Log_{10}( \rho [g][cm^{-3}] )$",fontsize=fs)
plt.legend()
plt.savefig("density_interpolated.png") #save
plt.show()

#Part 3: Solving ODE

#Get boundary conditions
mass_func = interp1d(r,mass,kind="cubic") #Interpolate mass data
R = r_grid[-1]
mass_R = mass_func(R) #Get M(r) at r=R_outer
phi_B = -ggrav*mass_R/R #Outer boundary condition
phi_A = phi(0,dens_grid[0],R) #Inner boundary condition

#3.2 Define functions for integration
def dudr(dens): return -4*np.pi*ggrav*dens

def calc_rhs(u,density): return u, dudr(density)

def integrate_FE(phi0,u0,density,dr,resolution):
    
    #Create array to return (phi,dphi) array
    yy = np.zeros( (resolution, 2) )
    
    #Use initial guesses
    yy[0,0] = phi0
    yy[0,1] = u0 
    print "%10s %10s %10s %10s %10s" % ("phi","phi'","dens","phi'","phi''")
    for i in range(resolution-1):
        #phi(n+1) = phi(n) + dr*phi'(n)
        #phi'(n+1) = phi'(n) + dr*phi''(n)
        u,du = calc_rhs(yy[i,1],density[i])
        print "%10.3e %10.3e %10.3e %10.3e %10.3e" % (yy[i,0],yy[i,1],density[i],u,du)
        yy[i+1,0] = yy[i,0] + dr*u
        yy[i+1,1] = yy[i,1] + dr*du
    
    return yy


#3.3 Create Variables and Initial Guesses
dr = r[1]-r[0] #step size
u0_old = 200.0 #Initial guess 1 for derivative of phi at r=0
u0_new = 10000.0 #Initial guess 2 for derivative of phi at r=0
yy0 = integrate_FE(phi_A,u0_old,dens_grid,dr,res) #profile 1
yy1 = integrate_FE(phi_A,u0_new,dens_grid,dr,res) #profile 2
phi0 = yy0[res-1,0] - phi_B #error in profile 1
phi1 = yy1[res-1,0] - phi_B #error in profile 2
dphidz = (phi1-phi0)/(u0_new-u0_old) #dphi/dz rate of change of 'z' wrt


i = 0 #Iteration counter
itmax = 100 #Max iterations
err = 1.0e99 #Error placeholder
criterion = 1.0e-17 #Desired accuracy
errs = [] #Array for errors
print u0_new,u0_old,phi1,phi0,dphidz
u0_old = u0_new #Shift newest value to 'old'
phi0 = phi1

while (err > criterion and i < itmax):
    print -phi0/dphidz
    u0_new = u0_old - phi0/dphidz #Secand update
    yy = integrate_FE(phi_A,u0_new,dens_grid,dr,res) #Get new array of values
    phi1 = yy[res-1,0] - phi_B #
    dphidz = (phi1-phi0)/(u0_new-u0_old) # dphi/dz numerical
    err = np.abs((u0_new-u0_old)/u0_old) # your error measure
    errs.append(err)
    print u0_new,u0_old,phi1,phi0,dphidz
    u0_old = u0_new
    phi0 = phi1
    i = i+1 


plt.figure()
plt.plot(r_grid,yy[:,0],label="Numerical")
plt.plot(r_grid,phi(r_grid,dens_grid,R),label="Analytical")
plt.legend(loc=4)
plt.xlabel(r"$r$",fontsize=fs)
plt.ylabel(r"$\phi(r)$",fontsize=fs)
plt.savefig("potential.png")
plt.show()


