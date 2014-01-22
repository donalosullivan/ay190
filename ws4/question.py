import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
#Function to find root of
def f(E,w,t,e): return E - w*t - e*np.sin(E)
#Analytical derivative
def df(E,e): return 1-e*np.cos(E)

def newton(old,w,t,e,count,f,df):
    new = old - f(old,w,t,e)/df(old,e)
    if abs(new-old) < 1e-10: return new,count
    else: return newton(new,w,t,e,count+1,f,df)
    
    
    
T = float(raw_input("Enter period (days): "))
a = float(raw_input("Enter semi-major axis (km):"))
e = float(raw_input("Enter eccentricity (e): "))
b = a*np.sqrt(1-e*82) #Calculate minor axis
w = 2*np.pi/T #angular frequency in rad/day

X = []
Y = []
print "%10s %10s %10s %10s %10s" % ("x","y","r","E","count")

for t in range(0,365):
    
    E_est = w*t/(1-e) #Small angle approximation for estimate/first guess
    E,count = newton(E_est,w,t,e,0,f,df)
    
    x = a*np.cos(E)
    y = b*np.sin(E)
    X.append(x)
    Y.append(y)
    r = np.sqrt( x**2 + y**2)
    print "%10.2f %10.2f %10.2f %10.2E %10i" % (x,y,r,E,count)
    
plt.figure()
plt.plot(X,Y,label="Orbit")
plt.legend()
plt.show()
