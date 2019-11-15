
import matplotlib.pyplot as plt
plt.switch_backend('agg')

# import matplotlib
# matplotlib.use('Agg')
import sys
sys.path.insert(0, "..")
import crnt4sbml
import numpy
import pandas
import sympy
import scipy.integrate as itg
from plotnine import ggplot, aes, geom_line, ylim, scale_color_distiller, facet_wrap, theme_bw, geom_path, geom_point
import time
import dill


network = crnt4sbml.CRNT("../sbml_files/GuanyuWang.xml")
# network.print_biological_reaction_types()
#
# ldt = network.get_low_deficiency_approach()
# ldt.report_deficiency_zero_theorem()
# ldt.report_deficiency_one_theorem()

#optimization approach
opt = network.get_mass_conservation_approach()
#opt.generate_report()


# the decision vector
opt.get_decision_vector()

# this function suggests physiological bounds
bounds, concentration_bounds = opt.get_optimization_bounds()

# overwriting specie concentration bounds for s4. Concentrations are in pM.
opt.get_concentration_bounds_species()

# params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds, iterations=1000,
#                                                                      concentration_bounds=concentration_bounds)


#numpy.save('params_guanyuwang.npy', params_for_global_min)

params_for_global_min = numpy.load('params_guanyuwang.npy')

# The reponse-related specie should be picked based on CellDesigner IDs. In our case phoshorylated A is s2.
# How to pick continuation parameter? In our case it is the amount of A protein, thus the conservation law 3.
print(opt.get_conservation_laws())

multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species="s7", parameters=params_for_global_min[[14]],
                                                           auto_parameters={'PrincipalContinuationParameter': 'C4'},
                                                           print_lbls_flag=True)

opt.generate_report()

sys.exit()


# odes = network.get_c_graph().get_ode_system()
# sympy_reactions = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_reactions()]
# sympy_species = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_species()]
# concentration_funs = opt.get_concentration_funs()
# BT_matrix = network.get_c_graph().get_b()
#
# important_variables = [odes, sympy_reactions, sympy_species, concentration_funs, BT_matrix, opt.get_decision_vector(),
#                        opt.get_conservation_laws(), multistable_param_ind, plot_specifications]
#
# dill.settings['recurse'] = True # allows us to pickle the lambdified functions
#
# with open("important_variables.dill", 'wb') as f:
#     dill.dump(important_variables, f)

with open("important_variables.dill", 'rb') as f:
    imprt_vars = dill.load(f)

odes = imprt_vars[0]
sympy_reactions = imprt_vars[1]
sympy_species = imprt_vars[2]
concentration_funs = imprt_vars[3]
BT_matrix = imprt_vars[4]
decision_vector = imprt_vars[5]
conservation_laws = imprt_vars[6]
multistable_param_ind = imprt_vars[7]
plot_specifications = imprt_vars[8]


# Parameters that produced bistability.
# re* are kinetic constants. Units can be found here help(network.get_physiological_range).
df = pandas.DataFrame(numpy.vstack([params_for_global_min[i] for i in multistable_param_ind]).T,
                      columns=["set" + str(i + 1) for i in multistable_param_ind],
                      index=[str(i) for i in decision_vector])
print(df)
#df.to_csv("succeeded_pars.txt", sep='\t')

################## selected parameter set #########################
#decision_vector_values = numpy.array(df['set4'])
# alternative declaration (for the sake of reference)
decision_vector_values = params_for_global_min[3]
plot_specifications = plot_specifications[0]  # warning, overwriting variable!!!

################ ODEs ###################################
print("Original ODEs")
#odes = network.get_c_graph().get_ode_system()
sympy.pprint(odes)

# String -> Sympy objects
# construct sympy form of reactions and species
#sympy_reactions = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_reactions()]
#sympy_species = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_species()]
# joining together
lambda_inputs = sympy_reactions + sympy_species
# creating a lambda function for each ODE to
ode_lambda_functions = [sympy.utilities.lambdify(lambda_inputs, odes[i]) for i in range(len(odes))]

############################### kinetic constants ########################################################
# Does this work for over, proper and under-dimensioned networks
kinetic_constants = numpy.array([decision_vector_values[i] for i in range(len(sympy_reactions))])

################################# Computing material conservation values ############################
# equilibrium species concentrations
species_concentrations = [i(*tuple(decision_vector_values)) for i in concentration_funs]
print(sympy_species)
print(species_concentrations)
print(conservation_laws)
# combine equilibrium specie concentrations according to conservation relationships
conservation_values = BT_matrix*sympy.Matrix([species_concentrations]).T

################################# starting concentrations ############################################
# this assumes that a chemical moiety in one state (specie) and other species containing this moiety are zero
# assignment of conservation values to species requires exploring the model in CellDesigner
# C1 is all in s4, free enzyme E2
# C2 is all in s3, free enzyme E1
# C3 is all in s1, free unphosphorylated specie A
# ['s2', 's3', 's1', 's1s2', 's7', 's7s2', 's4', 's5', 's7s4', 's5s3', 's6', 's9', 's9s6', 's8', 's3s8']
# [0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0]
# y_fwd = [conservation_values[4],  0.0,  conservation_values[3],  0.0,  0.0,                     0.0,  conservation_values[1],  0.0,  0.0,  0.0,  conservation_values[2],  0.0,  0.0,  conservation_values[0],  0.0]
# y_rev = [conservation_values[4],  0.0,  conservation_values[3],  0.0,  conservation_values[2],  0.0,  conservation_values[1],  0.0,  0.0,  0.0,  0.0,                     0.0,  0.0,  conservation_values[0],  0.0]

y_fwd = [conservation_values[4],  0.0,  conservation_values[3],  0.0,  0.0,    0.0,  conservation_values[1],  0.0,  0.0,  0.0,  conservation_values[2],  0.0,  0.0,  conservation_values[0],  0.0]
y_rev = [conservation_values[4],  0.0,  conservation_values[3],  0.0,  conservation_values[2],    0.0,  conservation_values[1],  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  conservation_values[0],  0.0]

print(y_fwd)
print(y_rev)

# Note, the continuation parameter C3 (first position) will be varied during simulations

############ simulation ###################
# computing dy/dt increments
# def f(cs, t, ks, ode_lambda_func, start_ind):
#     return [i(*tuple(ks), *tuple(cs)) for i in ode_lambda_func]  # dy/dt

def f(cs, t, ks, ode_lambda_func):
    return [i(*tuple(ks), *tuple(cs)) for i in ode_lambda_func]  # dy/dt

def sim_fun_fwd(x):
    y_fwd[2] = x
    return itg.odeint(f, y_fwd, t, args=(kinetic_constants, ode_lambda_functions))

def sim_fun_rev(x):
    y_rev[2] = x
    return itg.odeint(f, y_rev, t, args=(kinetic_constants, ode_lambda_functions))
#     return itg.odeint(f, y_fwd, t, args=(kinetic_constants, ode_lambda_functions, len(ode_lambda_functions)))21

# starting and ending time in seconds, number of data points
t = numpy.linspace(0.0, 3000000.0, 100000)
# signal parameter scanning range and data points. Forward scan.
# C3_scan = numpy.linspace(5.3e4, 5.4e4, 60)
# alternatively can be taken from plot_specifications
C3_scan = numpy.linspace(*plot_specifications[0], 30) #2.4e5,2.5e5, 30)  #
sim_res_fwd = [sim_fun_fwd(i) for i in C3_scan]  # occupies sys.getsizeof(sim_res_rev[0])*len(sim_res_rev)/2**20 Mb
# Reverse C3_scan. Reverse means that s2 is already high and signal is decreasing.
sim_res_rev = [sim_fun_rev(i) for i in numpy.flip(C3_scan)]

#one_C2 = [sim_res_rev[0][i][4] for i in range(len(sim_res_rev[0]))]
one_C2 = []
two_C2 = []

for i in range(len(C3_scan)):
    one_C2.append(sim_res_rev[i][-1][4])  #[4])
    two_C2.append(sim_res_fwd[i][-1][4]) #[4])

print(one_C2)
print(two_C2)

plt.plot(numpy.flip(C3_scan), one_C2, '.')
plt.plot(C3_scan, two_C2, '.')

#plt.plot(t, one_C2, '-')
#plt.ylim(0.0, 5000)
plt.savefig('./num_cont_graphs/bistability_ODEs.png')


print("hi")
################## exporting to text #####################################
# out = pandas.DataFrame(columns=['dir','signal','time'] + [str(i) for i in sympy_species])
# for i in range(len(sim_res_fwd)):
#     out_i = pandas.DataFrame(sim_res_fwd[i], columns=out.columns[3:])
#     out_i['time'] = t
#     out_i['signal'] = C3_scan[i]
#     out_i['dir'] = 'fwd'
#     out = pandas.concat([out, out_i[out.columns]])
# for i in range(len(sim_res_rev)):
#     out_i = pandas.DataFrame(sim_res_rev[i], columns=out.columns[3:])
#     out_i['time'] = t
#     out_i['signal'] = numpy.flip(C3_scan)[i]
#     out_i['dir'] = 'rev'
#     out = pandas.concat([out, out_i[out.columns]])
# out.to_csv("./num_cont_graphs/succeeded_pars.txt", sep="\t", index=False)
#
#
# ###################### plotting ##################################
# g = (ggplot(out, aes('time', 's7', group='signal', color='signal'))
#  + geom_line(size=0.5)
#  + ylim(0, plot_specifications[1][1])
#  + scale_color_distiller(palette='RdYlBu', type="diverging")
#  + facet_wrap('~dir')
#  + theme_bw())
# g.save(filename="./num_cont_graphs/sim_fwd_rev.png", format="png", width=8, height=4, units='in', verbose=False)
#
# eq = out[out.time == max(out.time)]
# g = (ggplot(eq)
#      + aes(x='signal', y='s7', color='dir')
#      + geom_path(size=2, alpha=0.5)
#      + geom_point(color="black")
#      + theme_bw())
# g.save(filename="./num_cont_graphs/sim_bif_diag.png", format="png", width=8, height=4, units='in', verbose=False)

