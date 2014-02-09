import numpy as np
import matplotlib.pyplot as plt
import time

#Gaussian Elimination
#Inputs: A - coefficient matrix, b - RHS vector of values
#Returns: x - solution vector
def gauss_elimination(A,b):

    n = A.shape[0] #Get dimensions
    if b.size!=n:
        print "Input dimensions do not work together"
        return np.zeros(n)
        
    #Sort so that pivots are non-zero: TBC
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
       
    #Set all 'close to zero' values to exactly 0.0       
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if abs(A[i][j]) < 1e-10:
                A[i][j]=0.0
    
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
    return M,BVEC
        
        
M,BVEC = get_data()
t1 = time.clock()
X = gauss_elimination(M[0],BVEC[0])
t2 = time.clock()
dt = t2-t1
print X,t1,t2,dt

