from scipy.optimize import minimize_scalar as minf
from matplotlib import pyplot as plt
import numpy as np

def f(x):
    return (1.0/4.0)*x**4.0 - (5.0/3.0)*x**3.0 - 6.0*x**2.0 + 19.0*x - 7
def fminus(x):
    return -(1.0/4.0)*x**4.0 + (5.0/3.0)*x**3.0 + 6.0*x**2.0 - 19.0*x + 7
def f2(x):
    return 3.0*x**4.0 - 4.0*x**3.0 + 1

#res = minf(fminus,bounds=(0,2),method='bounded')
#res = minf(f,bounds=(0,2),method='bounded')
res = minf(f2,bounds=(-4,3),method='bounded')
print res.x

'''
xx = np.linspace(-6,0,200)
plt.plot(xx,f(xx))
plt.show()
'''
