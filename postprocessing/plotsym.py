import sympy as sp
from sympy.plotting import plot3d
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc
plt.rc('text',usetex=True)

x, y = sp.symbols('x y')
#fun = x**2 + y**2
#fun = (1.0/2.0)*x**2.0 + (1.0/4.0)*y**4.0 - (1.0/2.0)*y**2.0
fun = (x-y**2)*(x-(1/2)*y**2)
numfun = sp.lambdify([x, y], fun)
xx, yy = np.linspace(-10,10,200), np.linspace(-10,10,200)
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
#ax.view_init(azim=0, elev=30)
plt.savefig('surf.eps')

plt.figure()
CS = plt.contour(X, Y, Zfun, 8, colors='k')
plt.clabel(CS, CS.levels, inline=True, inline_spacing=3, rightside_up=True, fontsize=8,fmt='%1.2f')
PC = plt.pcolor(X, Y, Zfun, cmap=cm.coolwarm)
clb = plt.colorbar(PC,shrink=0.5,aspect=5)
clb.ax.set_title(r'$f(x_1,x_2)$')
plt.xlabel(r'$x_1$')
plt.ylabel(r'$x_2$')
plt.savefig('contour.eps')
