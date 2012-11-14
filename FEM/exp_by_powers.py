import sympy as sm
from numpy import *
import matplotlib.pyplot as mpl

def least_squares(f,phi,Omega):
	N = len(phi)-1
	A = sm.zeros((N+1,N+1))
	b = sm.zeros((N+1,1))
	x = sm.Symbol('x')
	for i in range(N+1):
		for j in range(i,N+1):
			A[i,j] = sm.integrate(phi[i]*phi[j],(x,Omega[0],Omega[1]))
			A[j,i] = A[i,j]
		b[i,0] = sm.integrate(phi[i]*f,(x,Omega[0],Omega[1]))
	c = A.LUsolve(b)
	u = 0
	for i in range(len(phi)):
		u += c[i,0]*phi[i]
	return u

def comparison_plot(f,u,Omega,filename = 'tmp.pdf'):
	x = sm.Symbol('x')
	f = sm.lambdify([x],f,modules="numpy")
	u = sm.lambdify([x],u,modules="numpy")
	resolution = 401
	xcoor = linspace(Omega[0],Omega[1],resolution)
	exact = f(xcoor)
	approx = u(xcoor)
	mpl.plot(xcoor,approx)
	mpl.hold('on')
	mpl.plot(xcoor,exact)
	mpl.legend(['approximation','exact'])
	mpl.show()
	#savefig(filename)


x = sm.Symbol('x')
f = sm.exp(-x)
phi = [1,x,x**2,x**3,x**4,x**5,x**6,x**7,x**8]
Omega = [1,4]
u = least_squares(f,phi,Omega)
print u
comparison_plot(f,u,Omega)
