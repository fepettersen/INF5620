import sympy as sm
from numpy import *
import matplotlib.pyplot as mpl

def comparison_plot(f,u,Omega,filename = 'tmp.pdf'):
	x = sm.Symbol('x')
	f = sm.lambdify([x],f,modules="numpy")
	u = sm.lambdify([x],u,modules="numpy")
	resolution = 401
	#print f,phi,omega
	xcoor = linspace(Omega[0],Omega[1],resolution)
	exact = f(xcoor)
	approx = u(xcoor)
	mpl.plot(xcoor,approx)
	mpl.hold('on')
	mpl.plot(xcoor,exact)
	mpl.legend(['approximation','exact'])
	mpl.show()
	#savefig(filename)
	
def least_squares_orth(f,phi,Omega):
	N = len(phi)-1
	A = [0]*(N+1)
	b = [0]*(N+1)
	x = sm.Symbol('x')
	for i in range(N+1):
		A[i] = sm.integrate(phi[i]**2,(x,Omega[0],Omega[1]))
		b[i] = sm.integrate(phi[i]*f,(x,Omega[0],Omega[1]))
	c = [b[i]/A[i] for i in range(len(b))]
	u = 0
	print A,b
	for i in range(len(phi)):
		u+=c[i]*phi[i]
	return u

x = sm.Symbol('x')
phi = [sm.sin(x),sm.sin(2*x),sm.sin(3*x)]
omega = [0,3.14159*2]
f = sm.sin(20*x)

u = least_squares_orth(f,phi,omega)
print u
comparison_plot(f,u,omega)
