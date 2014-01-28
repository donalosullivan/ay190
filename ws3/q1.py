import numpy as np
from matplotlib import pyplot as plt

def f(x): return x*np.sin(x)

def int_f(a,b): return -np.cos(b)+np.cos(a)

def midpoint(a,b,fun,h):

    x = np.arange(a,b+h,h)
    integral = 0
    dx = x[1:] - x[:-1]
    for i in range(1,len(x)):
        # Qi = (bi - ai)*f(ai+bi/2)
        a,b = x[i-1],x[i]
        integral += ( b - a )*fun( (a+b)/2 )
        
    return integral
    
def trapezoid(a,b,fun,h):

    x = np.arange(a,b+h,h)
    integral = 0
    for i in range(1,len(x)):
        # Qi = (1/2)*(bi - ai)*{ f(a) + f(b) }
        a,b = x[i-1],x[i]
        integral += 0.5*( b - a )*(f(b) + f(a))
        
    return integral
    
def simpson(a,b,fun,h):
    x = np.arange(a,b+h,h)
    integral = 0
    for i in range(1,len(x)):
        # Qi = [ (b-a)/6 ]*{ f(a) + 4*f(a+b/2) + f(b) ]
        a,b = x[i-1],x[i]
        integral += (( b - a )/6)*(f(a) + 4*f((a+b)/2) + f(b))
        
    return integral
    
x1 = 0.0
x2 = np.pi
h = (x2-x1)/10
print "\nTesting convergence for numerical integration.",
print "h = %.1f, function = sin(x)" % h
print "\n%10s %10s %10s %10s" % ("","Midpoint","Trapezoid","Simpson")

hs = []
trap_errs = []
midp_errs = []
simp_errs = []

for i in range(0,5):
    hi = h/(2**i)
    actual = int_f(x1,x2)
    int_midp = midpoint(x1,x2,f,hi)
    int_trap = trapezoid(x1,x2,f,hi)
    int_simp = simpson(x1,x2,f,hi)
    err_midp = abs(int_midp - actual)
    err_trap = abs(int_trap - actual)
    err_simp = abs(int_simp - actual)
    
    hs.append(hi)
    midp_errs.append(err_midp)
    trap_errs.append(err_trap)
    simp_errs.append(err_simp)
    print "%10s %10.7E %10.7E %10.7E" % ("h/"+str(2**i),err_midp, err_trap, err_simp)
    
plt.figure()
fs = 24;
plt.plot(np.log10(hs),np.log10(midp_errs),label="Midpoint")
plt.plot(np.log10(hs),np.log10(trap_errs),label="Trapezoid")
plt.plot(np.log10(hs),np.log10(simp_errs),label="Simpson")
plt.xlabel("Log h",fontsize=fs)
plt.ylabel("Log error",fontsize=fs)
plt.legend(loc=4)
plt.show()
#plt.savefig("question1fig.png")




