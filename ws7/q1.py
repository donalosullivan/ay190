import numpy as np
import matplotlib.pyplot as plt
import random

#Random number generators
random.seed(1)
sr = random.SystemRandom()


def estimate_pi(Nvals,rand,plotstyle='ro',lab=""):
    #Variables needed
    r = 1.0
    pi_vals = []

    for N in Nvals: #Sample size loop

        count_in = 0.0 #Count hits inside radius
        
        for x in range(int(N)): #Random data loop
            x = r*rand()
            y = r*rand()
            dist = np.sqrt(x**2 + y**2)
            if dist <= r: count_in+=1
            
        pi = 4*count_in/N #Estimate pi
        pi_vals.append(pi)
        print "%10.2e %10.6f %10.6e" % (N,pi,abs(pi-np.pi))
        
    pi_vals = np.array(pi_vals)
    err_vals = abs(pi_vals - np.pi)/np.pi
    plt.plot(np.log10(Nvals),err_vals,plotstyle,label=lab)

fs = 20 #Fontsize for plot
plt.figure()
plt.xlabel(r"$Log_{10}(N)$",fontsize=fs)
plt.ylabel("% Error",fontsize=fs)
Nvals = np.linspace(1e3,5e4,100)
estimate_pi(Nvals,np.random.rand,'ro',"numpy.random")
estimate_pi(Nvals,sr.random,'bo',"random.SystemRandom")
plt.legend()
plt.savefig("q1fig.png")
plt.show()
