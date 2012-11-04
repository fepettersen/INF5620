# HPL: made loops over s, N and type of points.

import numpy as np
import sympy as sm
import matplotlib.pyplot as mpl

#u = lambda x: 4*L**L/3 - 4*L**2/np.pi*np.sin(np.pi*x/L)
'''
mpl.plot(x,f(x))
mpl.hold('on')
mpl.plot(x,u(x))
mpl.legend(['exact','approx'])
mpl.show()
'''
def Lagrange_polynomial(x, i, points):
	p = 1
	for k in range(len(points)):
		if k != i:
			p *= (x - points[k])/(points[i] - points[k])
	return p

def Lagrange_polynomials_01(x, N):
	if isinstance(x, sm.Symbol):
		h = sm.Rational(1, N-1)
	else:
		h = 1.0/(N-1)
	points = [i*h for i in range(N)]
	phi = [Lagrange_polynomial(x, i, points) for i in range(N)]
	return phi, points

def Chebyshev_nodes(a, b, N):
	from math import cos, pi
	return [0.5*(a+b) + 0.5*(b-a)*cos(float(2*i+1)/(2*(N+1))*pi) \
	for i in range(N+1)]

def Lagrange_polynomials(x, N, Omega, point_distribution='uniform'):
	if point_distribution == 'uniform':
		if isinstance(x, sm.Symbol):
			h = sm.Rational(Omega[1] - Omega[0], N)
		else:
			h = (Omega[1] - Omega[0])/float(N)
		points = [Omega[0] + i*h for i in range(N+1)]
	elif point_distribution == 'Chebyshev':
		points = Chebyshev_nodes(Omega[0], Omega[1], N)
	phi = [Lagrange_polynomial(x, i, points) for i in range(N+1)]
	return phi, points



def comparison_plot(f,u,Omega,filename = 'tmp.pdf', plot_title=''):
	x = sm.Symbol('x')
	f = sm.lambdify([x],f,modules="numpy")
	u = sm.lambdify([x],u,modules="numpy")
	resolution = 401
	#print f,phi,omega
	xcoor = np.linspace(Omega[0],Omega[1],resolution)
	exact = f(xcoor)
	approx = u(xcoor)
	mpl.plot(xcoor,approx)
	mpl.hold('on')
	mpl.plot(xcoor,exact)
	mpl.legend(['approximation','exact'])
	mpl.title(plot_title)
	mpl.show()
	mpl.savefig(filename)

def interpolation(f, phi, points):
	N = len(phi) - 1
	A = sm.zeros((N+1, N+1))
	b = sm.zeros((N+1, 1))
	x = sm.Symbol('x')
	# Turn phi and f into Python functions
	phi = [sm.lambdify([x], phi[i]) for i in range(N+1)]
	f = sm.lambdify([x], f)
	for i in range(N+1):
		for j in range(N+1):
			A[i,j] = phi[j](points[i])
		b[i,0] = f(points[i])
	c = A.LUsolve(b)
	u = 0
	for i in range(len(phi)):
		u += c[i,0]*phi[i](x)
	return u

x = sm.Symbol('x')
Omega = [0,1]
for s in 10,100:
    f = -sm.tanh(s*(x-0.5))
    for N in 3,6,9,11:
        for distribution in 'uniform', 'Chebyshev':
	    phi, points = Lagrange_polynomials(x, N, Omega,
					       point_distribution=distribution)

	    u = interpolation(f, phi, points)
	    filename = 'tmp_N_%d_s_%d_%s.pdf' % (N, s, distribution)
	    comparison_plot(f,u,Omega,filename,
			    plot_title='s=%g, N=%d, %s points' %
			    (s, N, distribution))
