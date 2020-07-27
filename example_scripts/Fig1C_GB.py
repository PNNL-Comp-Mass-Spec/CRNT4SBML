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

replacements = [[s1, (C3 - 2.0*s15 - s16 - s3 - s6)]]

ds2 = sympy.simplify(ds2.subs(s1, (C3 - 2.0*s15 - s16 - s3 - s6)))
ds3 = sympy.simplify(ds3.subs(s1, (C3 - 2.0*s15 - s16 - s3 - s6)))
ds6 = sympy.simplify(ds6.subs(s1, (C3 - 2.0*s15 - s16 - s3 - s6)))
ds15 = sympy.simplify(ds15.subs(s1, (C3 - 2.0*s15 - s16 - s3 - s6)))


# k1, k2, k3, k4, s, y, x, v0 = symbols('k1 k2 k3 k4 s y x v0', real=True)
# xdot = 2*k1*s*y - k2*x**2 - k3*x*y - k4*x + v0
# ydot = k2*x**2 - k1*s*y
# gb = groebner([xdot, ydot], [y, x], order='lex')
# print(gb)

# sys.exit()


# ds3 = re1*s1*s2 + s3*(-re1r - re2)
# ds6 = re2*s3 - re3*s6*s7 + re3r*s16 - re5*s1*s6 + s15*(re5r + 2*re6)
# ds16 = re3*s6*s7 + s16*(-re3r - re4)
# ds15 = re5*s1*s6 + s15*(-re5r - re6)

# ds15 = sympy.simplify(sympy.expand(re5*s6*(C3 - 2.0*s15 - s16 - s3 - s6) + s15*(-re5r - re6)))
# ds3 = sympy.simplify(sympy.expand(re1*(C2 - s3)*(C3 - 2.0*s15 - s16 - s3 - s6) + s3*(-re1r - re2)))
# ds6 = sympy.simplify(sympy.expand(re2*s3 - re3*s6*(C1 - s16) + re3r*s16 - re5*s6*(C3 - 2.0*s15 - s16 - s3 - s6) + s15*(re5r + 2.0*re6)))
# ds16 = sympy.simplify(sympy.expand(re3*s6*(C1 - s16) + s16*(-re3r - re4)))


# in CoCoLib:
# use QQab ::= QQ[re1, re1r, re2, re3, re3r, re4, re5, re5r, re6, C1, C2, C3];
# K := NewFractionField(QQab);
# use K[s3, s6, s16, s15], lex;
# I := ideal(C2*C3*re1 - 2.0*C2*re1*s15 - C2*re1*s16 - C2*re1*s3 - C2*re1*s6 - C3*re1*s3 + 2.0*re1*s15*s3 + re1*s16*s3 + re1*s3^2 + re1*s3*s6 - re1r*s3 - re2*s3, -C1*re3*s6 - C3*re5*s6 + re2*s3 + re3*s16*s6 + re3r*s16 + 2.0*re5*s15*s6 + re5*s16*s6 + re5*s3*s6 + re5*s6^2 + re5r*s15 + 2.0*re6*s15, C1*re3*s6 - re3*s16*s6 - re3r*s16 - re4*s16, C3*re5*s6 - 2.0*re5*s15*s6 - re5*s16*s6 - re5*s3*s6 - re5*s6^2 - re5r*s15 - re6*s15);
# SetVerbosityLevel(100);
# ReducedGBasis(I);

# SageMath Fig1Ci:
# A.<re1, re1r, re2, re3, re3r, re4, re5, re5r, re6, C1, C2, C3> = PolynomialRing(QQ)
# F = A.fraction_field()
# F.inject_variables()
# R.<s3, s6, s16, s15> = PolynomialRing(F, order='lex')
# I = R.ideal(C2*C3*re1 - 2.0*C2*re1*s15 - C2*re1*s16 - C2*re1*s3 - C2*re1*s6 - C3*re1*s3 + 2.0*re1*s15*s3 + re1*s16*s3 + re1*s3^2 + re1*s3*s6 - re1r*s3 - re2*s3, -C1*re3*s6 - C3*re5*s6 + re2*s3 + re3*s16*s6 + re3r*s16 + 2.0*re5*s15*s6 + re5*s16*s6 + re5*s3*s6 + re5*s6^2 + re5r*s15 + 2.0*re6*s15, C1*re3*s6 - re3*s16*s6 - re3r*s16 - re4*s16, C3*re5*s6 - 2.0*re5*s15*s6 - re5*s16*s6 - re5*s3*s6 - re5*s6^2 - re5r*s15 - re6*s15)
# I.groebner_basis(algorithm='singular:groebner', prot=True)

# Macaulay Fig1Ci:
# K = frac(QQ[re1, re1r, re2, re3, re3r, re4, re5, re5r, re6, C1, C2, C3]);
# R = K[s3, s6, s16, s15, MonomialOrder => Lex];
# gbTrace = 3
# J=ideal(C2*C3*re1 - 2*C2*re1*s15 - C2*re1*s16 - C2*re1*s3 - C2*re1*s6 - C3*re1*s3 + 2*re1*s15*s3 + re1*s16*s3 + re1*s3^2 + re1*s3*s6 - re1r*s3 - re2*s3, -C1*re3*s6 - C3*re5*s6 + re2*s3 + re3*s16*s6 + re3r*s16 + 2*re5*s15*s6 + re5*s16*s6 + re5*s3*s6 + re5*s6^2 + re5r*s15 + 2*re6*s15, C1*re3*s6 - re3*s16*s6 - re3r*s16 - re4*s16, C3*re5*s6 - 2*re5*s15*s6 - re5*s16*s6 - re5*s3*s6 - re5*s6^2 - re5r*s15 - re6*s15)
# G = ideal groebnerBasis(J, Strategy=>"F4")

# SageMath simple example:
# A.<k1, k2, k3, k4> = PolynomialRing(QQ)
# F = A.fraction_field()
# F.inject_variables()
# R.<x, y> = PolynomialRing(F, order='lex')
# I = R.ideal(2*k1*y-k2*x^2-k3*x*y-k4*x, k2*x^2-k1*y)
# I.groebner_basis(algorithm='singular:groebner', prot=True)

# Macaulay simple example:  MonomialOrder => Lex
# K = frac(QQ[k1, k2, k3, k4]);
# R = K[x,y, MonomialOrder => Lex];
# J=ideal(2*k1*y-k2*x^2-k3*x*y-k4*x, k2*x^2-k1*y)
# G = ideal groebnerBasis(J, Strategy=>"MGB");
# G = ideal groebnerBasis(J, Strategy=>"F4");
# netList G_*

# sys.exit()

def s_polynomial(f, g):
    return expand(lcm(LM(f), LM(g))*(1/LT(f)*f - 1/LT(g)*g))

from sympy.polys.polytools import parallel_poly_from_expr
import itertools


# species = (s1, s2, s3, s6, s7, s16, s15) # AA* is last
species = (s2, s3, s6, s7, s16, s15) # AA* is last
# species = [s3, s6, s16, s15]

ds2, ds3, ds6, ds7, ds15, ds16 = parallel_poly_from_expr([ds2, ds3, ds6, ds7, ds15, ds16], order='lex', gens=species,
                                               domain=RR[re1, re2, re3, re4, re5, re6, re1r, re3r, re5r, C3])[0]

# equations = [ds1, ds2, ds3, ds6, ds7, ds15, ds16]
equations = [ds2, ds3, ds6, ds7, ds15, ds16]

# print("Equations:")
# for i in equations:
#     print(i)
#     print("")
#
# sys.exit()

# equations = [ds3, ds6, ds16, ds15]

# equations_indices = [i for i in range(len(equations))]
# indices_combos = list(itertools.combinations(equations_indices, 2))
#
# all_spec_ind_combos = list(itertools.permutations(equations_indices, 4))
#
# all_spec_combos = []
# for i in all_spec_ind_combos:
#     all_spec_combos.append([species[ii] for ii in i])
#
# all_spec_combos = [(s3, s6, s16, s15), (s3, s16, s6, s15), (s6, s3, s16, s15), (s6, s16, s3, s15),
#                    (s16, s3, s6, s15), (s16, s6, s3, s15)]
#
# print(all_spec_combos)

# print("Starting")
# for ii in range(len(all_spec_combos)):
#     print("top of loop:")
#     print(f"index of species combos = {ii}")
#     ds3, ds6, ds16, ds15 = parallel_poly_from_expr([ds3, ds6, ds16, ds15], order='lex', gens=all_spec_combos[ii],
#                                                    domain=RR[s1, s2, s7, re1, re2, re3, re4, re5, re6,
#                                                              re1r, re3r, re5r, C1, C2, C3])[0]
#
#     F = [ds3, ds6, ds16, ds15]
#     iters = 0
#
#     while iters < 4:
#         F_indices = [i for i in range(len(F))]
#         indices_combos = list(itertools.combinations(F_indices, 2))
#         count = 0
#         for i in indices_combos:
#             s_poly = s_polynomial(F[i[0]], F[i[1]])
#             _, r = sympy.reduced(s_poly, F)
#             simp_expr = sympy.simplify(r.as_expr())
#             if simp_expr != sympy.S.Zero:
#                 F.append(sympy.simplify(r))
#                 break
#             else:
#                 count += 1
#         iters += 1
#
#         # if this conditions is met then a Groebner basis has already been obtained
#         if len(indices_combos) == count:
#             break
#
#     print(f"length of F = {len(F)}")
#     print(f"length of indices_combos = {len(indices_combos)}")
#     print(f"number of zero remainders = {count} \n")
#     print("\n \n")


print("---------------------------------------------------")

# #
import time

# ds3, ds6, ds16, ds15 = parallel_poly_from_expr([ds3, ds6, ds16, ds15], order='lex', gens=(s16, s3, s6, s15),
#                                                domain=RR[s1, s2, s7, re1, re2, re3, re4, re5, re6, re1r, re3r,
#                                                          re5r, C1, C2, C3])[0]
# # species = (s1, s2, s3, s6, s7, s16, s15) # AA* is last
# species = [s16, s3, s6, s15]
# # equations = [ds1, ds2, ds3, ds6, ds7, ds15, ds16]
# equations = [ds3, ds6, ds16, ds15]

start = time.time()
gb = groebner(equations, species, order='lex') #, method='f5b')
end = time.time()
elapsed = end - start
print(f"Elapsed time = {elapsed} \n")
print(gb)
print("Univariate basis polynomial = \n")
print(gb[-1])

# for i in range(len(gb)):
#     print(f"i={i}")
#     pprint(gb[i], num_columns=300)

# print("")
# univariate = gb[3]
# print("")
# replacements = [[s1, (C3 - 2.0*s15 - s16 - s3 - s6)], [s7, (C1 - s16)], [s2, (C2 - s3)]]
#
# for i in replacements:
#     univariate = univariate.subs(i[0], i[1])
#
# print("univariate with subs")
# expanded = sympy.expand(univariate)
# collected = sympy.collect(expanded, s15)
#
# pprint(collected)
