import sympy
from sympy import *
import sys
init_printing(num_columns=150, wrap_line=False, forecolor="red")

re1, re1r, re2, re3, re3r, re4, re5, re5r, re6 = symbols('re1 re1r re2 re3 re3r re4 re5 re5r re6', real=True)
s1, s2, s3, s6, s7, s15, s16 = symbols('s1 s2 s3 s6 s7 s15 s16', real=True)
C1, C2, C3 = symbols('C1 C2 C3', real=True)
ds1 = -re1*s1*s2 + re1r*s3 + re4*s16 - re5*s1*s6 + re5r*s15
ds2 = -re1*s1*s2 + s3*(re1r + re2)
ds3 = re1*s1*s2 + s3*(-re1r - re2)
ds6 = re2*s3 - re3*s6*s7 + re3r*s16 - re5*s1*s6 + s15*(re5r + 2*re6)
ds7 = -re3*s6*s7 + s16*(re3r + re4)
ds16 = re3*s6*s7 + s16*(-re3r - re4)
ds15 = re5*s1*s6 + s15*(-re5r - re6)




k1, k2, k3, k4, s, y, x, v0 = symbols('k1 k2 k3 k4 s y x v0', real=True)
xdot = 2*k1*s*y - k2*x**2 - k3*x*y - k4*x + v0
ydot = k2*x**2 - k1*s*y
gb = groebner([xdot, ydot], [y, x], order='lex')
print(gb)

sys.exit()




species = (s1, s2, s3, s6, s7, s16, s15) # AA* is last
# species = (s2, s3, s6, s7, s16, s15) # AA* is last

ds1, ds2, ds3, ds6, ds7, ds15, ds16 = parallel_poly_from_expr([ds1, ds2, ds3, ds6, ds7, ds15, ds16], order='lex', gens=species,
                                                              domain=RR[re1, re2, re3, re4, re5, re6, re1r, re3r, re5r])[0]

equations = [ds1, ds2, ds3, ds6, ds7, ds15, ds16]
# equations = [ds2, ds3, ds6, ds7, ds15, ds16]


import time

start = time.time()
gb = groebner(equations, species, order='lex', method='f5b')
end = time.time()
elapsed = end - start
print(f"Elapsed time = {elapsed} \n")
print(gb)
print("Univariate basis polynomial = \n")
print(gb[-1])

