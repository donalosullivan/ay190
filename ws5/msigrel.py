#!/usr/bin/env python

from astropy.io import ascii
from numpy import log10,arange,array
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline


def lin_reg(xdata,ydata,xerr=None,yerr=None):
    a1,a2 = -99,-99 #Random initial values for results
    Ex,Ey,Ex2,S,Exy = 0,0,0,0,0 #Declaring variables to use
    xerrdata = False if xerr==None else True
    yerrdata = False if yerr==None else True
    
    if xerrdata and yerrdata: 
        c,m = lin_reg(xdata,ydata) #Get slope of error-free data
        dx = xdata[1:]-xdata[:-1]
        dy = ydata[1:]-ydata[:-1]
        dydx = []
        for i in range(len(dx)):
            if dx[i]!=0: dydx.append(dy[i]/dx[i])
            else: dydx.append(m)
        dydx.append(dydx[-1])
        dydx = array(dydx) 
        yerr_ex = dydx*xerr
        sig_sqd = yerr**2 + yerr_ex**2 
    elif yerrdata: sig_sqd = yerr**2
    else: sig_sqd = array([1 for x in range(len(xdata))])
  
    Ex = sum(xdata/sig_sqd)
    Ey = sum(ydata/sig_sqd)
    Ex2 = sum((xdata**2)/sig_sqd)
    Exy = sum((xdata*ydata)/sig_sqd)
    S = sum(1/sig_sqd)
    a1 = (Ey*Ex2 - Ex*Exy)/(S*Ex2 - Ex**2)
    a2 = (S*Exy - Ey*Ex)/(S*Ex2 - Ex**2)
    print a1,a2
    return a1,a2
    
def line(xdata,a1,a2): return a1 + a2*xdata

#Open data file
data = ascii.read("m_sigma_table.dat",readme="m_sigma_ReadMe.dat")  

#Raw data is in form Log(M) and lin(sigma)
logM = array(data["logM"])
E_logM = array(data["E_logM"])
e_logM = array(data["e_logM"])
sig = array(data["sigma*"])
e_sig = array(data["e_sigma*"])
#Sort data for usability
inds = sig.argsort()
sig = sig[inds[::-1]]
logM = logM[inds[::-1]]
e_logM = e_logM[inds[::-1]]
e_sig = e_sig[inds[::-1]]

(sig_p,sig_m) = (sig + e_sig, sig - e_sig)
#Convert to log(sigma) for plotting purposes
log_sig = log10(sig)
log_sig_e = (log10(sig)-log10(sig_m),log10(sig_p)-log_sig)




#Non-error-based fit
a1,a2 = lin_reg(log_sig,logM)
fit = line(log_sig, a1, a2)  
fs = 18
plt.figure()
plt.subplot(211)
plt.errorbar(log_sig,
             logM,
             xerr=log_sig_e,
             yerr=e_logM,
             marker="o",
             color="k",
             linestyle="None")
plt.plot(log_sig,fit,label="No Error Fit: y=%.2fx + %.2f" % (a2,a1))
plt.ylabel(r'$Log_{10}M$',fontsize=fs)
plt.legend()
#Error-based fit
log_sig_e_avg = (abs(log_sig_e[0])+abs(log_sig_e[1]))/2
a1,a2 = lin_reg(log_sig,logM,log_sig_e_avg,e_logM) 
fit = line(log_sig, a1, a2) 

plt.subplot(212)
plt.errorbar(log_sig,
             logM,
             xerr=log_sig_e,
             yerr=e_logM,
             marker="o",
             color="k",
             linestyle="None")
plt.plot(log_sig,fit,label="Error Fit: y=%.2fx + %.2f" % (a2,a1))
plt.xlabel(r'$Log_{10}\sigma$',fontsize=fs)
plt.ylabel(r'$Log_{10}M$',fontsize=fs)
plt.legend()


plt.savefig("msigfit.png")
