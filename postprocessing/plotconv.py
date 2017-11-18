# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
import sys


plt.rc('text',usetex=True)
#plt.rc('font',family='helvica')
'''
###################################################
filename =  'mgrad1.dat'
data = np.genfromtxt(filename)
x, y = data[:,0], data[:,3]
#plt.xlim(xmin=1,xmax=20)
plt.plot(x,y,'-o',label='Grad')

###################################################
filename =  'mnewton1.dat'
data = np.genfromtxt(filename)
x, y = data[:,0], data[:,3]
#plt.xlim(xmin=1,xmax=20)
plt.plot(x,y,'-o',label='Newton')
'''
###################################################
filename =  'mdfp1.dat'
data = np.genfromtxt(filename)
x, y = data[:,0], data[:,3]
#plt.xlim(xmin=1,xmax=20)
plt.plot(x,y,'-o',label='DFP')

###################################################
filename =  'mdfp2.dat'
data = np.genfromtxt(filename)
x, y = data[:,0], data[:,3]
#plt.xlim(xmin=1,xmax=20)
plt.plot(x,y,'-o',label='sDFP')

###################################################
filename =  'mbfgs1.dat'
data = np.genfromtxt(filename)
x, y = data[:,0], data[:,3]
#plt.xlim(xmin=1,xmax=20)
plt.plot(x,y,'-o',label='BFGS')

###################################################
filename =  'mbfgs2.dat'
data = np.genfromtxt(filename)
x, y = data[:,0], data[:,3]
#plt.xlim(xmin=1,xmax=20)
plt.plot(x,y,'-o',label='sBFGS')

###################################################
filename =  'mdfp_bfgs1.dat'
data = np.genfromtxt(filename)
x, y = data[:,0], data[:,3]
#plt.xlim(xmin=1,xmax=20)
plt.plot(x,y,'-o',label='DFP/BFGS')

###################################################
filename =  'mdfp_bfgs2.dat'
data = np.genfromtxt(filename)
x, y = data[:,0], data[:,3]
#plt.xlim(xmin=1,xmax=20)
plt.plot(x,y,'-o',label='sDFP/sBFGS')
###################################################

plt.grid(True)
plt.xlabel(r'Iterations')
plt.ylabel(r'$f(x_1,x_2)$')
plt.legend(loc='best',borderpad=0.5)
plt.xlim(xmin=1,xmax=15)
plt.ylim(ymin=-0.1,ymax=10)
plt.savefig('plotconv.pdf')
