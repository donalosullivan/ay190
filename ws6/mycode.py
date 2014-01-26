from numpy import exp,array,pi
i = 1j
def dft(X):
    N = len(X)
    w = exp(-2*pi*i/N)
    A = array([ [w**(j*k) for j in range(0,N) ] for k in range(0,N) ])
    Y = A.dot(X)
    return Y
