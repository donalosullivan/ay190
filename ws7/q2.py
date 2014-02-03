from numpy.random import random_integers as randint
import numpy as np
import matplotlib.pyplot as plt

P_mcs,P_ans = [],[]
X=1000
Nvals = np.arange(1,30,1)

print "%10s %10s %10s" % ("N","P_MC","P_An")

for N in Nvals:

    #Monte carlo
    truecount = 0 #Count positive results
    
    for x in range(X): #Run 'X' simulations of size 'N'
    
        bdays = {} #Keep track of previous birthdays
        
        double = False
        
        for p in range(N): #Generate N birthdays 
            
            d = randint(1,365) 
            
            if bdays.has_key(d): #If already had, mark true and break
                double = True
                break
                
            else: #If new birthday, note and continue
                bdays[d] = True
                
        if double: #If this result produced a double, increment count
            truecount += 1
    P_mc = float(truecount)/X
    
    #Analytical
    P_an = 1
    for n in range(366-N,365): P_an *= n/365.0
    P_an = 1-P_an
    
    print "%10i %10.3f %10.3f" % (N,P_mc,P_an)
        
    P_mcs.append(P_mc)
    P_ans.append(P_an)


P_mcs,P_ans = np.array(P_mcs),np.array(P_ans)
print len(P_mcs),len(P_ans),len(Nvals)
   
fs = 20
plt.figure()
plt.xlabel("N",fontsize=fs)
plt.ylabel("P",fontsize=fs)
plt.plot(Nvals,P_mcs,'r-',label="Monte Carlo")
plt.plot(Nvals,P_ans,'b-',label="Analytic")
plt.plot(Nvals[22],P_mcs[22],'ko')
plt.text(Nvals[22]+2,P_mcs[22]-0.05,"P = %.2f\nN = %i" % (P_mcs[22],Nvals[22]))   
plt.legend()
plt.savefig("q2fig.png")
plt.show()      
        
    
