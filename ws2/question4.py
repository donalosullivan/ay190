import numpy as np
from matplotlib import pyplot as plt

#****************************
# LAGRANGE INTERPOLATION
#****************************

#Equation III.3.7
def L(n,j,x,x_):

    if j!=0: k = 0
    else: k = 1
    result = (x-x_[k])/(x_[j]-x_[k])
    for i in range(k,n):
        k += 1
        if k!=j:
            result *= (x-x_[k])/(x_[j]-x_[k])   
    return result
        
#Equation III.3.6
def p(n,x,x_,y_):

    result = 0
    for j in range(0,n):
        result += y_[j]*L(n,j,x,x_)
    return result

#****************************
# PIECEWISE INTERPOLATION
#****************************

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

def getinterval2(x,x_):
    closest = 1e10
    a,b = 0,1
    for i,xi in enumerate(x_):
        if xi==x: return -1,i,0 #special case, x is a data point
        if abs(x-xi)<closest:
            closest = abs(x-xi)
            if x<xi: 
                a=i-1
                b=i
            else:
                a=i
                b=i+1
    c = 0
    if a==0: c=2 #if in first interval, use first 3 as abc
    elif b==len(x_)-1: #if in last interval, use last 3
        c == len(x_)-1
        b -= 1
        a -= 1
    else:
        c == b-1 #If not touching either boundary, use a,b,c
    
    return a,b,c
    
#Linear - Equation III.3.3
def pw_lin(x,x_,y_):
    a,b = getinterval(x,x_)
    if a==-1: return y_[b] #special case, x is a data point
    return y_[a] + ( (y_[b] - y_[a])/(x_[b] - x_[a]) )*(x - x_[a])

#Quadratic - Equation III.3.4
def pw_quad(x,x_,y_):
    a,b,c = getinterval2(x,x_)
    if a==-1: return y_[b] #special case, x is a data point
    x0,x1,x2 = x_[a],x_[b],x_[c]
    f0,f1,f2 = y_[a],y_[b],y_[c]
    result = (((x-x1)*(x-x2))/((x0-x1)*(x0-x2)))*f0
    result += (((x-x0)*(x-x2))/((x1-x0)*(x1-x2)))*f1
    result += (((x-x0)*(x-x1))/((x2-x0)*(x2-x1)))*f2
    return result
    
time = np.array([0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0])
appmag = np.array([0.302,0.185,0.106,0.093,0.24,0.579,0.561,0.468,0.302])

x = np.arange(0.0,1.0,0.001)
p_ = p(8,x,time,appmag)
p_pwlin = []
for xi in x: p_pwlin.append(pw_lin(xi,time,appmag))
p_pwquad = []
for xi in x: p_pwquad.append(pw_quad(xi,time,appmag))

fs=24
plt.figure()
plt.plot(x,p_,'g-',label="n=8 Lagrange")
plt.xlabel("time (days)",fontsize=fs)
plt.ylabel("apparent magnitude",fontsize=fs)

plt.plot(x,p_pwlin,'r-',label="Piecewise Linear")
plt.plot(x,p_pwquad,'b-',label="Piecewise Quadratic")
plt.plot(time,appmag,'go',label="data")
plt.legend(loc=9)
plt.savefig("question4fig.png")




