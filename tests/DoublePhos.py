import matplotlib.pyplot as plt
#plt.switch_backend('agg')

# import matplotlib
# matplotlib.use('Agg')

import crnt4sbml
import numpy
import pandas
import sympy
import scipy.integrate as itg
from plotnine import ggplot, aes, geom_line, ylim, scale_color_distiller, facet_wrap, theme_bw, geom_path, geom_point
import time

# import sys
# sys.path.insert(0, "..")

network = crnt4sbml.CRNT("../sbml_files/DoublePhos.xml")

# network.print_biological_reaction_types()
#
# ldt = network.get_low_deficiency_approach()
# ldt.report_deficiency_zero_theorem()
# ldt.report_deficiency_one_theorem()

# optimization approach
opt = network.get_mass_conservation_approach()
# opt.generate_report()


# the decision vector
#opt.get_decision_vector()

# this function suggests physiological bounds
bounds, concentration_bounds = opt.get_optimization_bounds()

# overwriting specie concentration bounds for s4. Concentrations are in pM.
opt.get_concentration_bounds_species()

params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds,
                                                                     concentration_bounds=concentration_bounds)

# The reponse-related specie should be picked based on CellDesigner IDs. In our case phoshorylated A is s2.
# How to pick continuation parameter? In our case it is the amount of A protein, thus the conservation law 3.
#print(opt.get_conservation_laws())

multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species="s4", parameters=params_for_global_min,
                                                           auto_parameters={'PrincipalContinuationParameter': 'C2'},
                                                           print_lbls_flag=False)

opt.generate_report()



#
# numpy.save('./Double_Phos/the_params.npy', params_for_global_min)
#
# import dill
#
# odes = network.get_c_graph().get_ode_system()
# sympy_reactions = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_reactions()]
# sympy_species = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_species()]
# concentration_funs = opt.get_concentration_funs()
# BT_matrix = network.get_c_graph().get_b()
#
# important_variables = [odes, sympy_reactions, sympy_species, concentration_funs, BT_matrix, multistable_param_ind,
#                        plot_specifications]
#
# dill.settings['recurse'] = True # allows us to pickle the lambdified functions
#
# with open("./Double_Phos/important_variables.dill", 'wb') as f:
#     dill.dump(important_variables, f)

# import dill
#
# with open("./Double_Phos/important_variables.dill", 'rb') as f:
#     out_1 = dill.load(f)
#
# params_for_global_min = numpy.load('./Double_Phos/the_params.npy')
#
# odes = out_1[0]
#
# # construct sympy form of reactions and species
# sympy_reactions = out_1[1]
# sympy_species = out_1[2]
# plot_specifications = out_1[6]


# Parameters that produced bistability.
# re* are kinetic constants. Units can be found here help(network.get_physiological_range).
df = pandas.DataFrame(numpy.vstack([params_for_global_min[i] for i in multistable_param_ind]).T,
                      columns=["set" + str(i + 1) for i in multistable_param_ind],
                      index=[str(i) for i in opt.get_decision_vector()])
print(df)

################## selected parameter set #########################
#decision_vector_values = numpy.array(df['set3'])
# alternative declaration (for the sake of reference)
decision_vector_values = params_for_global_min[5]
plot_specifications = plot_specifications[5]  # warning, overwriting variable!!!

################ ODEs ###################################
print("Original ODEs")
odes = network.get_c_graph().get_ode_system()
sympy.pprint(odes)

# String -> Sympy objects
# construct sympy form of reactions and species
sympy_reactions = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_reactions()]
sympy_species = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_species()]
# joining together
lambda_inputs = sympy_reactions + sympy_species
# creating a lambda function for each ODE to
ode_lambda_functions = [sympy.utilities.lambdify(lambda_inputs, odes[i]) for i in range(len(odes))]

############################### kinetic constants ########################################################
# Does this work for over, proper and under-dimensioned networks
kinetic_constants = numpy.array([decision_vector_values[i] for i in range(len(network.get_c_graph().get_reactions()))])

################################# Computing material conservation values ############################
# equilibrium species concentrations
species_concentrations = [i(*tuple(decision_vector_values)) for i in opt.get_concentration_funs()]
print(network.get_c_graph().get_species())
print(species_concentrations)
print(opt.get_conservation_laws())
# combine equilibrium specie concentrations according to conservation relationships
conservation_values = network.get_c_graph().get_b()*sympy.Matrix([species_concentrations]).T

################################# starting concentrations ############################################
# this assumes that a chemical moiety in one state (specie) and other species containing this moiety are zero
# assignment of conservation values to species requires exploring the model in CellDesigner
# C1 is in s4, free enzyme E2
# C2 is in s3, free enzyme E1
# C3 is in s1, free unphosphorylated specie A
# ['s1', 's2', 's3', 's3s1', 's4', 's4s2', 's2s1']
# ['C3',    0, 'C2',      0, 'C1',      0,      0]


# ['s1', 's5', 's2', 's2s1', 's3', 's3s5', 's4', 's2s5', 's3s4']
# ['C3',  0,   'C2',    0,   'C1',     0,     0,     0,   0]

y_fwd = [conservation_values[2], 0.0, conservation_values[1], 0.0, conservation_values[0], 0.0, 0.0, 0.0, 0.0]


y_rev = [conservation_values[2], 0.0, conservation_values[1], 0.0, conservation_values[0], 0.0, 0.0, 0.0, 0.0]

# y_fwd = [conservation_values[2], 0.0, conservation_values[1], 0.0, conservation_values[0], 0.0, 0.0]
# y_rev = [0.0, conservation_values[2], conservation_values[1], 0.0, conservation_values[0], 0.0, 0.0]
# Note, the continuation parameter C3 (first position) will be varied during simulations



############ simulation ###################
# computing dy/dt increments
def f(cs, t, ks, ode_lambda_func, start_ind):
    return [i(*tuple(ks), *tuple(cs)) for i in ode_lambda_func]  # dy/dt

def sim_fun_fwd(x):
    y_fwd[0] = x  # updating s1 concentration or C3
    return itg.odeint(f, y_fwd, t, args=(kinetic_constants, ode_lambda_functions, len(ode_lambda_functions)))

def sim_fun_rev(x):
    y_rev[1] = x  # updating s2 concentration
    return itg.odeint(f, y_rev, t, args=(kinetic_constants, ode_lambda_functions, len(sympy_reactions)))

# starting and ending time in seconds, number of data points
t = numpy.linspace(0.0, 300.0, 300)
# signal parameter scanning range and data points. Forward scan.
# C3_scan = numpy.linspace(5.3e4, 5.4e4, 60)
# alternatively can be taken from plot_specifications
C3_scan = numpy.linspace(*plot_specifications[0], 30)
print(C3_scan)
sim_res_fwd = [sim_fun_fwd(i) for i in C3_scan]  # occupies sys.getsizeof(sim_res_rev[0])*len(sim_res_rev)/2**20 Mb

#print("hi")
# Reverse C3_scan. Reverse means that s2 is already high and signal is decreasing.
sim_res_rev = [sim_fun_rev(i) for i in numpy.flip(C3_scan)]

#print("hi hi")
################## exporting to text #####################################
out = pandas.DataFrame(columns=['dir','signal','time'] + network.get_c_graph().get_species())
for i in range(len(sim_res_fwd)):
    out_i = pandas.DataFrame(sim_res_fwd[i], columns=out.columns[3:])
    out_i['time'] = t
    out_i['signal'] = C3_scan[i]
    out_i['dir'] = 'fwd'
    out = pandas.concat([out, out_i[out.columns]])
for i in range(len(sim_res_rev)):
    out_i = pandas.DataFrame(sim_res_rev[i], columns=out.columns[3:])
    out_i['time'] = t
    out_i['signal'] = numpy.flip(C3_scan)[i]
    out_i['dir'] = 'rev'
    out = pandas.concat([out, out_i[out.columns]])

#print("yo")
out.to_csv("doublephos_sim.txt", sep="\t", index=False)

#print("yo yo")
###################### plotting ##################################
g = (ggplot(out, aes('time', 's4', group='signal', color='signal'))
 + geom_line(size=0.5)
 + ylim(0, plot_specifications[1][1])
 + scale_color_distiller(palette='RdYlBu', type="diverging")
 + facet_wrap('~dir')
 + theme_bw())
#print("yo yo yo")
g.save(filename="./num_cont_graphs/sim_fwd_rev.png", format="png", width=8, height=4, units='in', verbose=False)
#print("ya ya")
eq = out[out.time == max(out.time)]
g = (ggplot(eq)
     + aes(x='signal', y='s4', color='dir')
     + geom_path(size=2, alpha=0.5)
     + geom_point(color="black")
     + theme_bw())
#print("ya")
g.save(filename="./num_cont_graphs/sim_bif_diag.png", format="png", width=8, height=4, units='in', verbose=False)
#print("ya yo")
