"""
Refresh linalg skills.
Chapter: Finite element methods, Exercise 1
"""

#Show that vectors in the plane (a,b) from a vectorspace by showing that all 
#axioms of a vector space is satisfied

#systematic solution:
#1: u+v = v+u
u+v = (a,b)+(c,d) = (a+c,b+d) = (c+a,b+d) = v+u #because a,b,c,d \in R
u+v = ax+b + cx+d = (a+c)x +b+d = (c+a)x +d+b = v+u #for linear functions

#2: u+(v+w) = (u+v)+w
u+(v+w) = (a,b) + (c+e,d+f) = (a+c+e,b+d+f) = (a+c,b+d) + (e,f) = (u+v)+w #same reason as in 1
u+(v+w) = ax+b + ((c+e)x +d+f) = (a+c+e)x +b+d+f = ((a+c)x +b+d) + ex+f = (u+v)+w # for linear functions

#3: u+0 = u
u+0 = (a,b) + (0,0) = (a+0,b+0) = (a,b) = u
u+0 = ax+b + 0x+0 = (a+0)x + b + 0 = ax+b = u

#4: u-u = 0
u-u = (a,b)-(a,b) = (a-a,b-b) = (0,0) = 0
u-u = ax+b - (ax+b) = (a-a)x +b -b = 0x+0 = 0

#5: k(u+v) = ku+kv
k(u+v) = k(a+c,b+d) = (k(a+c),k(b+d)) = (ka+kc,kb+kd) = (ka,kc) + (kc,kd) = ku+kv
k(u+v) = k((a+c)x +b+d) = k(a+c)x +kb+kd = kax+kb + kcx + kd = ku + kv

#6: (a+b)v = av + bv
(a+b)v = ((a+b)c,(a+b)d) = (ac+bc,ad+bd) = av + bv
(a+b)v = (a+b)cx + (a+b)d = acx +ad + bcx + bd = av + bv

#7: a(bv) = (ab)v 
a(bv) = a(bc,bd) = (abc,abd) = (ab)v
a(bv) = a(bcx + bd) = abcx + abd = (ab)v

#8 1*v = v
1*v = (1*c,1*d) = (c,d) = v
1*v = 1*cx + 1*d = cx+d = v

"""
Chapter: Finite element methods, Exercise 2
"""


