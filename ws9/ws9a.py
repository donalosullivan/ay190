import numpy as np
import matplotlib.pyplot as plt
import time
from numpy.linalg import solve as npsolve
from scipy.linalg import solve as spsolve

#Gaussian Elimination
#Inputs: A - coefficient matrix, b - RHS vector of values
#Returns: x - solution vector
def gauss_elimination(A,b):

    n = A.shape[0] #Get dimensions
    if b.size!=n:
        print "Input dimensions do not work together"
        return np.zeros(n)
        
    #Sort so that pivots are non-zero: Ignored for now
    for i in range(n):
        for j in range(n):
            pass

    #Put in upper triangular form
    for k in range(n-1):
    
        for j in range(k+1,n):
            
            ajk = float(A[j][k])
            akk = float(A[k][k])
            A[j] = -(ajk/akk)*A[k] + A[j]
            b[j] = -(ajk/akk)*b[k] + b[j]
       
    
    #Back substitution to solve for x values
    X = np.zeros(n)
    for j0 in range(n):
        j = n-1-j0
        X[j] = b[j]
        for i in range(j+1,n): X[j] -= A[j][i]*X[i]  
        X[j] /= A[j][j]
        
    return X
            

def get_data():
    M,BVEC = [],[]

    print "Loading data..."

    for i in range(1,6):   
        m = np.loadtxt("LSE%i_m.dat" % i)
        bvec = np.loadtxt("LSE%i_bvec.dat" % i)
        M.append(m)
        BVEC.append(bvec)

    print "Data loaded. Checking for solvability"

    solvable = [True for i in range(5)]
    for i in range(len(M)):
        
        det = np.linalg.slogdet(M[i])[1]
        if det==0.0: solvable[i]=False

    for i in range(len(solvable)):
        if not solvable[i]:
            print "Cannot solve %i" % i
        else:
        
            m = M[i]
            bvec = BVEC[i]
    return M,BVEC,solvable

def timer(solver,A,b):
    #Use time.time because time.clock doesn't work with Ubuntu
    t1 = time.time()
    X = solver(A,b)
    t2 = time.time()
    return t2-t1

def part1():
    
    A,b,solvable = get_data()
    for i in range(len(b)):
        print "N=",A[i].shape[0],
        if solvable[i]: print "Solvable"
        else: print "Unsolvable"
        
def part2():
    A = np.array([[2.0,3.0,4.0],[1.0,-1.0,1.0],[1.0,2.0,-1.0]])
    b = np.array([29.5,7.5,-3.0])
    #Known solution: [x1,x2,x3] = [2.0, 0.5, 6.0]
    X = gauss_elimination(A,b)
    X2 = npsolve(A,b)
    print "My solver: ",X
    print "NumPy solver: ",X2
    
def part3(): 
    A,b,solvable = get_data() #Pull data from files
    N = [] #Array to store x-axis values of problem size
    np_linalg_dt = [] #Numpy Linalg Solve times
    sp_linalg_dt = [] #Scipy Linalg Solve times
    my_linalg_dt = [] #Gauss_elimination solve times
    for i in range(len(b)): #Run through data
        print i+1,"/",len(b)
        if not solvable[i]: continue
        N.append(A[i].shape[0]) 
        np_linalg_dt.append(timer(npsolve,A[i],b[i])) #Time Numpy
        sp_linalg_dt.append(timer(spsolve,A[i],b[i])) #Time Scipy
        my_linalg_dt.append(timer(gauss_elimination,A[i],b[i])) #Time my solver
    N = np.array(N)
    #Log values
    np_linalg_dt = np.log10(np_linalg_dt)
    sp_linalg_dt = np.log10(sp_linalg_dt)
    my_linalg_dt = np.log10(my_linalg_dt)
    #Plot
    fs=20
    plt.figure()
    plt.plot(N,np_linalg_dt,'o',label="Numpy Linalg")
    plt.plot(N,sp_linalg_dt,'o',label="Scipy Linalg")
    plt.plot(N,my_linalg_dt,'o',label="Gauss Elim")
    plt.legend(loc=4)
    plt.xlabel(r"$N$",fontsize=fs)
    plt.ylabel(r"$\Delta t$",fontsize=fs)
    plt.savefig("fig1.png")
    plt.show()    
    
part1() #Print sizes of arrays
part2() #Test solver on known system
#part3() #Compare with NumPy and SciPy solvers
