import sys
sys.path.insert(0, "..")
import crnt4sbml
import numpy

# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/a_b.xml") # yes 10
# signal = "C1"
# response = "s5"
iters = 50

network = crnt4sbml.CRNT("../sbml_files/Fig1Ci.xml")
signal = "C3"
response = "s15"

network.basic_report()
network.print_c_graph()

GA = network.get_general_approach()
GA.initialize_general_approach(signal=signal, response=response, fix_reactions=True)


import sympy
# sympy.pprint(network.get_c_graph().get_s())
# print("")
# sympy.pprint(network.get_c_graph().get_s().rref())
# sys.exit()

print(GA.get_conservation_laws())
print(GA.get_fixed_reactions())
print(network.get_c_graph().get_reactions())
import sympy
from sympy import *
print(network.get_c_graph().get_ode_system())
print(network.get_c_graph().get_species())
# species = GA.get_independent_species()
# print(GA.get_independent_species())

# re5, re5r, re6, re6r, re7, re7r = sympy.symbols('re5, re5r, re6, re6r, re7, re7r', positive=True)
#
# # sympy.pprint(network.get_c_graph().get_a())
# #
print("")
print(GA.get_independent_odes_subs())
species = GA.get_independent_species()
print(GA.get_independent_species())
# indp_subs = GA.get_independent_odes_subs()

import itertools

s1, s2, s3, s6, s7, s15, s16 = symbols('s1 s2 s3 s6 s7 s15 s16', real=True)

species = [s1, s2, s3, s6, s7, s16, s15]
all_spec_combos = list(itertools.permutations(species, 7))

# print(all_spec_combos)
for i in all_spec_combos:
    if i[6] == s15:
        print(i)

sys.exit()
#
# C1 = sympy.symbols('C1', real=True)
#
# ds5, ds7 = sympy.parallel_poly_from_expr([indp_subs[0], indp_subs[1]], order='lex', gens=species,
#                                           domain=sympy.RR[re5, re5r, re6, re6r, re7, re7r, C1])[0]


# print(ds5)
# print("")
# print(ds7)


print("")


# import itertools
#
#
# def buchberger(F):
#     """Toy implementation of Buchberger algorithm. """
#
#     def s_polynomial(f, g):
#         return expand(lcm(LM(f), LM(g)) * (1 / LT(f) * f - 1 / LT(g) * g))
#
#     remainders = F
#     while sum(remainders) != sympy.S.Zero:
#         remainders = []
#         combos = list(itertools.combinations(range(len(F)), 2))
#
#         for i in combos:
#
#
#
#
#
#     return None
#
# F = [indp_subs[0], indp_subs[1]]
# buchberger(F)

# equations = [ds5, ds7]
#
# gb = sympy.groebner(equations, species, order='lex', method='f5b')
#
# print(gb[-1])

#
# sympy.pprint(GA.get_solutions_to_fixed_reactions())
#
# sys.exit()
# print(network.get_c_graph().get_species())
# print(network.get_c_graph().get_ode_system())
# print("")
# print(GA.get_independent_species())
# print(GA.get_independent_odes())
# print("")
# print(GA.get_independent_odes_subs())

bnds = [(1e-2, 1e2)]*len(GA.get_input_vector())
# params_for_global_min, obj_fun_vals = GA.run_optimization(bounds=bnds, iterations=iters, seed=0, print_flag=False,
#                                                           dual_annealing_iters=1000, confidence_level_flag=True, parallel_flag=True)
#
# if GA.get_my_rank() == 0:
#     numpy.save('./basic_example_paper/params.npy', params_for_global_min)
params_for_global_min = numpy.load('./basic_example_paper/params.npy')

# print(GA.get_input_vector())
# print(params_for_global_min[1])

multistable_param_ind, plot_specifications = GA.run_greedy_continuity_analysis(species=response, parameters=params_for_global_min,
                                                                               auto_parameters={'PrincipalContinuationParameter': signal},
                                                                               dir_path="./basic_example_paper", plot_labels=["$B_{tot}$", "[A]", None])

# list_of_ggplots = GA.run_direct_simulation([params_for_global_min[1]], parallel_flag=True, print_flag=True)
#
# if GA.get_my_rank() == 0:
#
#     print("list_of_ggplots ")
#     g = list_of_ggplots[0]
#     import plotnine as p9
#     path = "./basic_example_paper"
#     from matplotlib import rc
#     rc('text', usetex=True)
#
#     g = (g + p9.xlab("$B_{tot}$") + p9.ylab("$[A]$") + p9.scale_color_hue(labels=["High [A]", "Low [A]"]))
#     g.save(filename=path + f"/sim_bif_diag_0_0.png", format="png", width=6, height=4, units='in', verbose=False)
#
#     print("")

# GA.get_ode_lambda_functions()
#
# GA.get_jac_lambda_function()


def ff(t, cs, ks, ode_lambda_functions, jacobian):
    return [i(*tuple(ks), *tuple(cs)) for i in ode_lambda_functions]

def jac_f(t, cs, ks, ode_lambda_functions, jacobian):
    return jacobian(*tuple(ks), *tuple(cs))


# forward scan:
# s5 = 0.0
# s6 = 0.0
# s7 = C1
#
# reverse scan :
# s5 = 0.0
# s6 = C1
# s7 = 0.0

param_index = 1

print(params_for_global_min[param_index])
sys.exit()

C1_val = params_for_global_min[param_index][7] + params_for_global_min[param_index][8]

y_fwd = [10.0, 0.0, C1_val]
y_rev = [0.0, C1_val, 0.0]

import scipy.integrate as itg
import pandas
from plotnine import ggplot, aes, geom_line, ylim, scale_color_distiller, facet_wrap, theme_bw, geom_path, geom_point, labs, annotate
from matplotlib import rc
rc('text', usetex=True)

kinetic_constants = params_for_global_min[param_index][0:6]
end_t = 65000.0
t = numpy.linspace(0.0, end_t, 32500)
scan_vals = numpy.linspace(7.77, 7.805, 60)

# end_t = 30000.0
# t = numpy.linspace(0.0, end_t, 100000)
# scan_vals = numpy.linspace(1.77, 1.790, 100)

def sim_fun_fwd(x):
    y_fwd[2] = x
    return itg.solve_ivp(ff, (0.0, end_t), y_fwd, t_eval=t, method="BDF", args=(kinetic_constants, GA.get_ode_lambda_functions(),
                                                       GA.get_jac_lambda_function()), jac=jac_f)

def sim_fun_rev(x):
    y_rev[1] = x
    return itg.solve_ivp(ff, (0.0, end_t), y_rev, t_eval=t, method="BDF", args=(kinetic_constants, GA.get_ode_lambda_functions(),
                                                       GA.get_jac_lambda_function()), jac=jac_f)

sim_res_fwd = [numpy.transpose(sim_fun_fwd(i).y) for i in scan_vals]
sim_res_rev = [numpy.transpose(sim_fun_rev(i).y) for i in numpy.flip(scan_vals)]

out = pandas.DataFrame(columns=['dir', 'signal', 'time'] + network.get_c_graph().get_species())
for i in range(len(sim_res_fwd)):
    out_i = pandas.DataFrame(sim_res_fwd[i], columns=out.columns[3:])
    out_i['time'] = t
    out_i['signal'] = scan_vals[i]
    out_i['dir'] = 'High [A]'
    out = pandas.concat([out, out_i[out.columns]])
for i in range(len(sim_res_rev)):
    out_i = pandas.DataFrame(sim_res_rev[i], columns=out.columns[3:])
    out_i['time'] = t
    out_i['signal'] = numpy.flip(scan_vals)[i]
    out_i['dir'] = 'Low [A]'
    out = pandas.concat([out, out_i[out.columns]])
# out.to_csv("./basic_example_paper/sim.txt", sep="\t", index=False)

###################### plotting ##################################
g = (ggplot(out, aes('time', response, group='signal', color='signal'))
     + geom_line(size=0.5)
     + ylim(0, 5)
     + labs(x="time", y="[A]")
     + scale_color_distiller(palette='RdYlBu', type="diverging", name="$B_{tot}$")
     + facet_wrap('~dir')
     + theme_bw())
g.save(filename="./basic_example_paper/sim_fwd_rev.png", format="png", width=8, height=4, units='in', verbose=False)

eq = out[out.time == max(out.time)]
g = (ggplot(eq)
     + aes(x='signal', y=response, color='dir')
     + labs(x="$B_{tot}$", y="[A]", color="")
     + geom_path(size=2, alpha=0.5)
     + geom_point(color="black")
     + theme_bw()
     + geom_point(color="black"))

g.save(filename="./basic_example_paper/sim_bif_diag.png", format="png", width=6, height=4, units='in', verbose=False)