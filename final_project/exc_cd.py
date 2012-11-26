from dolfin import *
import numpy as np

def alpha(u):
	return 1+u**2
maxiter = 20
rho = 3.14
dt = 0.05
mesh = UnitSquare(50,50)
V = FunctionSpace(mesh,'Lagrange',1)
f = Constant(0.0)
u = TrialFunction(V)
v = TestFunction(V)

u0 = Expression('cos(pi*x[0])',pi=pi)
u_1 = interpolate(u0,V)

def u0_boundary(x,on_boundary):
	return on_boundary

bc = DirichletBC(V,u0,u0_boundary)

a = u*v*dx + dt*inner(nabla_grad(u),nabla_grad(v))*dx
L = (u_1 +dt*f)*v*dx

A = assemble(a)
u = Function(V)
T = 2
t = dt
b = None

while t<=T:
	b = assemble(L, tensor=b)
	u0.t = t
	bc.apply(A,b)
	solve(A,u.vector(),b)
	plot(u)
	#rescale=False
	#interactive()
	t+=dt
	u_1.assign(u)
	u_e = interpolate(u0, V)
	maxdiff = np.abs(u_e.vector().array()-u.vector().array()).max()
	#print 'Max error, t=%.2f: %-10.17f' % (t, maxdiff)
