import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

def getinterval(x,x_):
    closest = 1e10
    a,b = 0,1
    for i,xi in enumerate(x_):
    
        if xi==x: return -1,i #special case, x is a data point
        if abs(x-xi)<closest:
            closest = abs(x-xi)
            if x<xi: 
                a=i-1
                b=i
            else:
                a=i
                b=i+1
    return a,b
    
#**************************************
# PIECEWISE CUBIC HERMITE INTERPOLATION
#**************************************

def df(x,y):
    result = []
    for i in range(1,len(x)):
        result.append( (y[i]-y[i-1])/(x[i]-x[i-1]))
    return result
        
        
def pw3Hermite(x,x_,y_):
    d = df(x_,y_)
    a,b = getinterval(x,x_)
    if a==-1: return y_[b]
    x1,x2 = x_[a],x_[b]
    y1,y2 = y_[a],y_[b]
    
    if a<0: return
    d1,d2 = d[a-1],d[b-1]
    
    z = (x-x1)/(x2-x1)

    def psi0(z): return 2*z**3 - 3*z**2 + 1
    def psi1(z): return z**3 - 2*z**2 + z
    
    H3 = y1*psi0(z) + y2*psi0(1-z) + d1*(x2-x1)*psi1(z) - d2*(x2-x1)*psi1(1-z)
    
    return H3
    
    

time = np.array([0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0])
appmag = np.array([0.302,0.185,0.106,0.093,0.24,0.579,0.561,0.468,0.302])

x = np.arange(0.0,1.01,0.01)
p_pw3H = []
for xi in x: p_pw3H.append(pw3Hermite(xi,time,appmag))

p_ncs = interp1d(time,appmag,kind='cubic')
fs=24
plt.figure();
plt.plot(x,p_pw3H,'r-',label="PW Cubic Hermite")
plt.plot(x,p_ncs(x),'b-',label="Natural Cubic Spline")
plt.xlabel("time (days)",fontsize=fs)
plt.ylabel("apparent magnitude",fontsize=fs)

plt.plot(time,appmag,'go',label="data")
plt.legend(loc=4)
plt.savefig("question5fig.png")

