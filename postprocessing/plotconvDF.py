# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
import sys


plt.rc('text',usetex=True)
#plt.rc('font',family='helvica')

###################################################
filename =  'mbox1.dat'
data = np.genfromtxt(filename)
x, y = data[:,0], data[:,3]
#plt.xlim(xmin=1,xmax=20)
plt.plot(x,y,'-o',label='Box')

###################################################
filename =  'mHJ1.dat'
data = np.genfromtxt(filename)
x, y = data[:,0], data[:,3]
#plt.xlim(xmin=1,xmax=20)
plt.plot(x,y,'-o',label='Hooke-Jeeves')

###################################################
plt.grid(True)
plt.xlabel(r'Iterations')
plt.ylabel(r'$f(x_1,x_2)$')
plt.legend(loc='best',borderpad=0.5)
plt.xlim(xmax=8)
#plt.ylim(ymin=-0.1)
plt.savefig('plotconv.eps')
