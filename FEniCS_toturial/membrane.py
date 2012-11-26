# coding=UTF-8
from dolfin import *
import numpy

T = 10.0
A = 1.0
R = 0.3
theta = 0.2
x0 = 0.6*R*cos(theta)
y0 = 0.6*R*sin(theta)
sigma = 50

mesh = UnitCircle(30)
V = FunctionSpace(mesh,'Lagrange',1)

#Define boundary conditions
u0 = Constant(0.0)
def u0_boundary(x,on_boundary):
	return on_boundary

bc = DirichletBC(V,u0,u0_boundary)

#define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression('4*exp(-0.5*pow((R*x[0]-x0)/sigma,2) -0.5*pow((R*x[1]-y0)/sigma,2))',x0=x0,R=R,y0=y0,sigma=sigma)
a = inner(nabla_grad(u),nabla_grad(v))*dx
L = f*v*dx

#Compute solution
u = Function(V)
problem = LinearVariationalProblem(a,L,u,bc)
solver = LinearVariationalSolver(problem)
solver.parameters['linear_solver'] = 'cg'
solver.parameters['preconditioner'] = 'ilu'
solver.solve()
'''
solve(a == L, u, bc,
	solver_parameters={'linear_solver': 'cg',
						'preconditioner': 'ilu'})
'''
#Plot solution and mesh
plot(u, wireframe=True, title='solution')
plot(mesh)

f = interpolate(f,V)
plot(f,title='scaled pressure')
max_u = u.vector().array().max()
max_D = A*max_u/(8*pi*sigma*T)

#Verification for "flat" pressure (large sigma)
if sigma>=50:
	u_e = Expression('1-x[0]*x[0] - x[1]*x[1]')
	u_e = interpolate(u_e,V)
	dev = numpy.abs(u_e.vector().array() - u.vector().array()).max()
	print 'sigma=%g: max deviation=%e' % (sigma, dev)

'''
file = File('poisson.pvd')
file << u
'''
#Hold plot
interactive()