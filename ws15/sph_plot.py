import numpy as np
import os
import matplotlib.pyplot as mpl

files = os.listdir(".")
dat_files = {}
for x in files:
    if ".dat" in x and "out_0." in x:
        try: num = int(x[6:8])
        except:
            try: num = int(x[6:7])
            except: continue
        dat_files[num] = x
        
mpl.figure()

for key,value in dat_files.items():
    print key,value
    i,r,rho,press,v = np.loadtxt(value).transpose()
    if key==0:
        mpl.subplot(221)
        mpl.plot(r,rho,'+',label="t = 0.%02i" % key)
        mpl.legend()
        mpl.ylabel(r"$\rho$",fontsize=20)
    elif key==6:
        mpl.subplot(222)
        mpl.plot(r,rho,'+',label="t = 0.%02i" % key)
        mpl.legend()
    elif key==12:
        mpl.subplot(223)
        mpl.plot(r,rho,'+',label="t = 0.%02i" % key)
        mpl.legend()
        mpl.ylabel(r"$\rho$",fontsize=20)
        mpl.xlabel(r"$r$",fontsize=20)
    elif key==18:
        mpl.subplot(224)
        mpl.plot(r,rho,'+',label="t = 0.%02i" % key)
        mpl.legend() 
        mpl.xlabel(r"$r$",fontsize=20)
mpl.savefig("time_plot.png")
mpl.show()
        
