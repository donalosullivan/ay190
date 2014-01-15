from matplotlib import pyplot as plt
import numpy as np

onethird = np.float32(1.0/3.0)

def x(n):

	if n==0: return np.float32(1)
	elif n==1: return onethird
	else: return np.float32(13*onethird*x(n-1) - 4*onethird*x(n-2))

def x_analytic(n): return onethird**n

print "Absolute error values:"
print "%10s"%"n", "%8s"%"Abs_Err","%8s"%"Rel_Err"


X = []
absolute_err = []
relative_err = []
for n in range(2,16):
	abs_err = abs(x(n)-x_analytic(n))
	rel_err = abs_err/x_analytic(n)
	X.append(n)
	absolute_err.append(np.log10(abs_err))
	relative_err.append(np.log10(rel_err))
	print "%10i"%n,"%6.3E" % abs_err,"%6.3E" % rel_err

plt.figure()
plt.title("Error in recurrence relation")
plt.plot(X,absolute_err,label="Abolute")
plt.plot(X,relative_err,label="Relative")
plt.ylabel("Log (Error magnitude)")
plt.xlabel("n")
plt.legend()
plt.savefig("error_plot.png")
plt.show()
