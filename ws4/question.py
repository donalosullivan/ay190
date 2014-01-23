import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

#Function to find root of
def f(E,w,t,e): return E - w*t - e*np.sin(E)

#Analytical derivative
def df(E,e): return 1-e*np.cos(E)

#Newton Raphson root finding
def newton(old,w,t,e,f,df,count):
    if count>100: #To avoid standard error output, handle iterations error manually
        print "Exceeded max iterations"
        return old,count
    if f(old,w,t,e)==0.0: return old,count #If exact solution, i.e f(x)=0, return
    change = - f(old,w,t,e)/df(old,e) #Get change in estimate
    new = old + change #Update guess
    fn0,fn1 = f(old,w,t,e),f(new,w,t,e)    
    frac_err = abs(( fn1 - fn0 )/fn0 ) #Get fraction error
    if frac_err < 1e-10: return new,count #If within our desired accuracy, return
    else: return newton(new,w,t,e,f,df,count+1) #Otherwise, next guess...
  
T = 365.25635 #Period in days 
a = 1.496e6 #Semi major axis in km 
e = float(raw_input("Eccentricity: ")) #Get eccentricity from user
b = a*np.sqrt(1-e**2) #Calculate minor axis
w = 2*np.pi/T #angular frequency in rad/day

X = []
Y = []

output = open("question1_%s.txt" % str(e), 'w')

output.write("%10s %10s %10s %10s %10s %10s %10s\n" % ("day","w*t","x (km)","y (km)","r (km)","E","count"))
print "%10s %10s %10s %10s %10s %10s %10s\n" % ("day","w*t","x (km)","y (km)","r (km)","E","count")


for t in [91,182,273]:
    
    E_est = w*t/(1-e) #Small angle approximation for estimate/first guess
    E,count = newton(E_est,w,t,e,f,df,0) #Get value for E by Newton-Raphson
    x = a*np.cos(E)
    y = b*np.sin(E)
    X.append(x)
    Y.append(y)
    r = np.sqrt( x**2 + y**2)
    output.write("%10i %10.2f %10.3e %10.3e %10.3e %10.3E %10i\n" % (t,w*t,x,y,r,E,count))
    print "%10i %10.2f %10.3e %10.3e %10.3e %10.3E %10i\n" % (t,w*t,x,y,r,E,count)
