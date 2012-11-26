from dolfin import *
import numpy as np
import sys

degree = int(sys.argv[1])
divisions = [int(i) for i in sys.argv[2:]]
d = len(divisions) -1
domain_type = [UnitInterval,UnitSquare,UnitCube]
mesh = domain_type[d](*divisions)

def alpha(u):
	return 1+u**2

def picard(u,u_1,a,L,bc,maxiter,tol=1e-5,order=2):
	iter = 0;eps = 1.0
	while(eps>tol and iter<maxiter):
		iter +=1
		solve(a==L,u,bc)
		diff = u.vector().array()-u_1.vector().array()
		eps = np.linalg.norm(diff,ord=order)
		u_1.assign(u)
	return None

maxiter = 20
rho = 3.14
dt = 0.1
V = FunctionSpace(mesh,'Lagrange',degree)
f = Constant(0.0)
u = TrialFunction(V)
v = TestFunction(V)

u0 = Expression('cos(pi*x[0])',pi=pi)
u_1 = interpolate(u0,V)

def u0_boundary(x,on_boundary):
	return on_boundary

bc = DirichletBC(V,u0,u0_boundary)

a = u*v*dx + dt*inner(alpha(u_1)*nabla_grad(u),nabla_grad(v))*dx
L = (u_1 +dt*f)*v*dx

A = assemble(a)
u = Function(V)
T = 1
t = dt
b = None
exact = Expression('exp(-pi*pi*t)*cos(pi*x[0])',pi=pi,t=0)
while t<=T:
	b = assemble(L, tensor=b)
	u0.t = t
	bc.apply(A,b)
	picard(u,u_1,a,L,bc,maxiter)
	#solve(A,u.vector(),b)
	plot(u)
	#rescale=False
	interactive()
	t+=dt
	u_1.assign(u)
	exact.t=t
	u_e = interpolate(exact, V)
	maxdiff = np.abs(u_e.vector().array()-u.vector().array()).max()
	#print 'Max error, t=%.2f: %-10.17f' % (t, maxdiff)
	e = u_e.vector().array() - u.vector().array()
	E = np.sqrt(np.sum(e**2)/u.vector().array().size)
