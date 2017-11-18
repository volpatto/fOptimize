import sympy as sp
from sympy.plotting import plot3d
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc
import sys
plt.rc('text',usetex=True)

'''
filename = sys.argv[1]
data = np.genfromtxt(filename)
xsol, ysol = data[:,1], data[:,2]
'''

x, y = sp.symbols('x y')
#fun = x**2 + y**2
#fun = (1.0/2.0)*x**2.0 + (1.0/4.0)*y**4.0 - (1.0/2.0)*y**2.0
#fun = (x-y**2)*(x-(1/2)*y**2)
#fun = (x**2 + y - 11)**2 + (x + y**2 -7)**2 # Himmelblau
fun = 100*(y-x**2)**2 + (x-1)**2 # Rosenbrock
numfun = sp.lambdify([x, y], fun)
inf = -5.0
sup = 5.0
npoints = 200
xx, yy = np.linspace(inf,sup,npoints), np.linspace(inf,sup,npoints)
X, Y = np.meshgrid(xx,yy)
Zfun = numfun(X,Y)

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Zfun, cmap=cm.coolwarm)
ax.set_xlabel(r'$x_1$')
ax.set_ylabel(r'$x_2$')
ax.set_zlabel(r'$f(x_1,x_2)$')
clb = fig.colorbar(surf,shrink=0.5,aspect=5)
clb.ax.set_title(r'$f(x_1,x_2)$')
#ax.view_init(azim=-70, elev=10)
plt.savefig('surf.pdf')

plt.figure()
CS = plt.contour(X, Y, Zfun, 15, colors='k')
plt.clabel(CS, CS.levels, inline=True, inline_spacing=3, rightside_up=True, fontsize=8,fmt='%1.2f')
PC = plt.pcolor(X, Y, Zfun, cmap=cm.coolwarm)
clb = plt.colorbar(PC,shrink=0.5,aspect=5)
clb.ax.set_title(r'$f(x_1,x_2)$')
plt.xlabel(r'$x_1$')
plt.ylabel(r'$x_2$')

'''
plt.plot(xsol, ysol, 'k', zorder=1, lw=1.5)
plt.scatter(xsol, ysol, s=20, color='k', zorder=2)
plt.xlim((inf,sup))
plt.ylim((inf,sup))
'''
plt.axes().set_aspect('equal')
#plt.plot(xsol,ysol,"-.",color='k')
plt.savefig('contour.pdf')
