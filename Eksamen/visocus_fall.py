#!/usr/bin/env python

from math import pi
import matplotlib.pyplot as mpl
from numpy import linspace,zeros
import sys, odespy

if len(sys.argv)!=4:
	print "usage: %s theta dt N"%sys.argv[0]

g = 9.81
rho = 1000
rho_f = 1.23
mu = 0.3
r = 0.01
m = (4.0/3)*pi*r**3*rho

theta = float(sys.argv[1])
dt = float(sys.argv[2])
N = int(sys.argv[3])

a = 6*pi*r*mu/m
b = rho_f/rho -1
v = [0]
t = [0]

def step_theta(v,a,b,dt,theta):
	return (v*(1-a*dt*(1-theta)) + b*dt)/(1+a*dt*theta)

def f(y,t):
	return b*g -a*y

t_d = linspace(0,N*dt,N+1)
solver = odespy.RK4(f)
solver.set_initial_condition(0.0)
v_rk,t_d = solver.solve(t_d)

for i in range(1,N):
	new = step_theta(v[i-1],a,b,dt,theta)
	v.append(new)
	t.append(i*dt)

#mpl.plot(t_d,v_rk)

mpl.plot(t,v)
mpl.ylabel('velocity')
mpl.xlabel('time in seconds')
mpl.show()