from dolfin import *

mesh = UnitSquare(5,5)
V_g = VectorFunctionSpace(mesh,'Lagrange',1)

w = TrialFunction(V_g)
v = TestFunction(V_g)

a = inner(w,v)*dx
L = inner(nabla_grad(w),nabla_grad(v))*dx

grad_u = Function(V_g)
solve(a==L,grad_u)

plot(grad_u,title='grad(u)')