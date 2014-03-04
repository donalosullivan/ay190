#!/usr/bin/env python

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mpl3d

# global constants
ggrav = 6.67e-8
msun  = 1.99e33
seconds_per_year = 24.*3600*365 # roughly
cm_per_pc = 3.1e18
distance_to_sgrAstar = 8e3 * cm_per_pc

# system parameters - UNCOMMENT to change file and units to sgrAstar.dat
initial_data_file = "sgrAstar.asc"
distance_unit_to_cm = 0.04*cm_per_pc
time_unit_to_s = seconds_per_year
mass_unit_to_g = msun
Nsteps = int(1e3)
t0 = 0
t1 = 100*seconds_per_year
dt = (t1-t0)/Nsteps

final_data_file = "final_positions.asc"

def NbodyRHS(u,mass):

    x,y,z = u[:,0],u[:,1],u[:,2]
    vx,vy,vz = u[:,3],u[:,4],u[:,5]
    
    n = len(x)
    
    #initialize empty accelerations
    ax,ay,az = np.zeros(n),np.zeros(n),np.zeros(n) 
    for i in range(n):
         
        for j in range(n): #sum through other particles and add grav accelerations   
              
            if i!=j:
                rij = np.array((x[i]-x[j],y[i]-y[j],z[i]-z[j])) #displacements
                Rij = np.sqrt(rij[0]**2 + rij[1]**2 + rij[2]**2)
                
                aij = (-ggrav*mass[j]/Rij**3)*rij 
                ax[i] += aij[0]
                ay[i] += aij[1]
                az[i] += aij[2]
                
                
    return np.array((vx,vy,vz,ax,ay,az)).transpose()    
                

def NbodyRK4(u,mass,dt):
    
    k1 = dt*NbodyRHS(u,mass)
    k2 = dt*NbodyRHS(u+k1/2,mass)
    k3 = dt*NbodyRHS(u+k2/2,mass)
    k4 = dt*NbodyRHS(u+k3,mass)

    return  u + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    
def TotalEnergy(u,mass):

    x,y,z = u[:,0],u[:,1],u[:,2]
    vx,vy,vz = u[:,3],u[:,4],u[:,5]
    
    v_squared = vx**2 + vy**2 + vz**2
    Ek = sum( 0.5*mass*v_squared) #kinetic energy
    
    n = len(x)
    Ep = 0 #potential energy
    for i in range(n):
        for j in range(n):
            if j!=i:
                rij = np.sqrt( (x[i]-x[j])**2 + (y[i]-y[j])**2 + (z[i]+z[j])**2 )
                Ep +=ggrav*mass[i]*mass[j]/rij
                
    return Ek + Ep
    

# main program
plt.ion()

(x,y,z,vx,vy,vz,mass) = np.loadtxt(initial_data_file, unpack = True)

# convert from units in initial data file to cgs
x *= distance_unit_to_cm
y *= distance_unit_to_cm
z *= distance_unit_to_cm
vx *= distance_unit_to_cm / time_unit_to_s
vy *= distance_unit_to_cm / time_unit_to_s
vz *= distance_unit_to_cm / time_unit_to_s
mass *= mass_unit_to_g

xmin = np.amin(x)
xmax = np.amax(x)
ymin = np.amin(y)
ymax = np.amax(y)
zmin = np.amin(z)
zmax = np.amax(z)
rmax = 2.5*max(abs(xmin),abs(xmax),abs(ymin),abs(ymax),abs(zmin),abs(zmax))

# use a single state vector to simplify the ODE code
# indices:
# u[:,0] = x
# u[:,1] = y
# u[:,2] = z
# u[:,3] = vx
# u[:,4] = vy
# u[:,5] = vz
u = np.array((x,y,z,vx,vy,vz)).transpose()

x,y,z,vx,vy,vz = u[0,:]
plt.clf()
fig = plt.gcf()
ax = mpl3d.Axes3D(fig)
ax.scatter(u[:,0],u[:,1],u[:,2])
ax.set_xlim((-rmax,rmax))
ax.set_ylim((-rmax,rmax))
ax.set_zlim((-rmax,rmax))
plt.draw()

energies,times = [],[]

for it in range(0, Nsteps):
    time = t0 + it * dt
    u = NbodyRK4(u,mass,dt)
    if it % max(1,Nsteps/10) == 0:
        E = TotalEnergy(u,mass)
        energies.append(E)
        times.append(time)
        print "\n"
        print "it = %d, time = %g years, energy = %g" % (it, time / seconds_per_year,E)
        plt.clf()
        fig = plt.gcf()
        ax = mpl3d.Axes3D(fig)
        ax.scatter(u[:,0],u[:,1],u[:,2])
        ax.set_xlim((-rmax,rmax))
        ax.set_ylim((-rmax,rmax))
        ax.set_zlim((-rmax,rmax))
        plt.draw()
       



# output result
file_header = "1:x 2:y 3:z 4:vx 5:vy 6:vz 7:mass"
np.savetxt(final_data_file, u, header=file_header)

energies,times = np.array(energies),np.array(times)
energies = np.log10(np.abs(energies - energies[0])/energies)

plt.figure()
plt.plot(times,energies)
plt.ylabel(r"$Log_{10} ( \frac{\Delta E}{E})$",fontsize=20)
plt.xlabel(r"$Time [s]$",fontsize=20)
plt.savefig("energy_sgrA.png")
plt.show()

