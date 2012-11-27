from dolfin import *
import numpy as np
import sys

degree = int(sys.argv[1])
divisions = [int(i) for i in sys.argv[2:]]
d = len(divisions) -1
domain_type = [UnitInterval,UnitSquare,UnitCube]
mesh = domain_type[d](*divisions)

def alpha(u):
	beta = 1.3
	return 1.0+beta*u**2

def picard(u,u_1,a,L,b,maxiter,tol=1e-5,order=2):
	iter = 0;eps = 1.0
	while(eps>tol and iter<maxiter):
		iter +=1
		#bc.apply(A,b)
		solve(A,u.vector(),b)
		diff = u.vector().array()-u_1.vector().array()
		eps = np.linalg.norm(diff,ord=order)
		u_1.assign(u)
	#print iter
	return None

maxiter = 1
rho = 1.0
dt = 0.01
sigma = 0.1

V = FunctionSpace(mesh,'Lagrange',degree)
#f = Expression('rho*x[0]*x[0]*(1.0/2-x[0]/3.0) +pow(x[0],4)*pow(t,3)*(8.0*pow(x[0],3)/9.0 -28.0*x[0]*x[0]/9.0 +\
#				7.0*x[0]/2.0 -5.0/4.0) +t*(2.0*x[0] -1.0)',rho=rho,t=0.0)#Constant(0.0)
f = Constant(0.0)
u = TrialFunction(V)
v = TestFunction(V)

u0 = Expression('exp(-1/(2*sigma*sigma)*(x[0]*x[0]+x[1]*x[1]))',sigma=sigma)
#u0 = Expression('cos(pi*x[0])',pi=pi)
#u0 = Constant(0.0)
u_1 = interpolate(u0,V)

def u0_boundary(x,on_boundary):
	return on_boundary

#bc = DirichletBC(V,u0,u0_boundary)

a = u*v*dx + dt*inner(alpha(u_1)*nabla_grad(u),nabla_grad(v))*dx
L = (u_1 +dt*f)*v*dx

A = assemble(a)
u = Function(V)
T = 0.5
t = dt
b = None
exact = Expression('t*x[0]*x[0]*(1.0/2- x[0]/3.0)',t=0.0)
plt_lst = [0.05,0.1,0.25,0.4]
counter=0
while t<=T:
	b = assemble(L, tensor=b)
	u0.t = t
	f.t=t
	picard(u,u_1,a,L,b,maxiter)
	#solve(A,u.vector(),b)
	
	#rescale=False
	#interactive()
	t+=dt
	u_1.assign(u)
	u_e = interpolate(exact, V)
	if counter ==25:
		viz_v = plot(u_e,title='exact',basename ='exact')
		viz_u = plot(u,title='nummeric',basename='nummeric')
		interactive()
		#viz_u.update(u)
		#viz_u.write_ps('nummeric',format='pdf')

	maxdiff = np.abs(u_e.vector().array()-u.vector().array()).max()
	exact.t=t
	#print 'Max error, t=%.2f: %-10.17f' % (t, maxdiff)
	e = u_e.vector().array() - u.vector().array()
	E = np.sqrt(np.sum(e**2)/u.vector().array().size)
	counter +=1
	if t==0.05:
		print E/dt

'''
exp(-pi*pi*t)*cos(pi*x[0])
'''