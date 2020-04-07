import sys
sys.path.insert(0, "..")
import crnt4sbml
import numpy
import sympy
import pandas
import scipy.integrate as itg
import dill
from plotnine import ggplot, aes, geom_line, ylim, scale_color_distiller, facet_wrap, theme_bw, geom_path, geom_point, labs, annotate
from matplotlib import rc
rc('text', usetex=True)

network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts_submodel_4c.xml")
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts_submodel_1.xml")
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts.xml")

signal = "C1"
response = "s11"

# network.basic_report()
# network.print_c_graph()

# GA = network.get_general_approach(signal=signal, response=response, fix_reactions=False)

GA = network.get_general_approach()
print(GA.get_conservation_laws())

GA.initialize_general_approach(signal=signal, response=response)

# print(GA.get_conservation_laws())

# submodel 4c bounds
# bnds = [(16.5, 17.5), (92.0, 92.5), (0.01, 0.025), (0.2, 0.25), (0.78, 0.79), (3.6, 3.7), (0.15, 0.25), (0.06, 0.065),    (0.9, 1.1), (1.9, 2.1)] + \
#        [(25.85, 25.95), (8.1, 8.15), (18.6, 18.67), (88.85, 88.95), (99.7, 99.8), (30.6, 30.7), (0.0, 1.0)]

bnds = [(16.5, 17.5), (92.0, 92.5), (0.01, 0.025), (0.2, 0.25), (0.78, 0.79), (3.6, 3.7), (0.15, 0.25), (0.06, 0.065),    (0.0, 100.0), (0.0, 100.0)] + \
       [(25.85, 25.95), (8.1, 8.15), (18.6, 18.67), (88.85, 88.95), (99.7, 99.8), (30.6, 30.7), (0.0, 100.0)]

# Nuts bounds
# bnds = [(1e-2, 100.0), (1e-2, 100.0), (1e-2, 100.0), (1e-2, 100.0), (0.9, 1.1), (1.9, 2.1), (16.5, 17.5), (92.0, 92.5), (0.01, 0.025), (0.2, 0.25), (0.78, 0.79), (3.6, 3.7), (0.15, 0.25), (0.06, 0.065)] + \
#        [(0.0, 1.0), (18.6, 18.67), (1e-2, 100.0), (1e-2, 100.0), (25.85, 25.95), (8.1, 8.15), (88.85, 88.95), (99.7, 99.8), (30.6, 30.7)]

# bnds = [(1e-2, 100.0), (1e-2, 100.0), (1e-2, 100.0), (1e-2, 100.0), (0.0, 100.0), (0.0, 100.0), (16.5, 17.5), (92.0, 92.5), (0.01, 0.025), (0.2, 0.25), (0.78, 0.79), (3.6, 3.7), (0.15, 0.25), (0.06, 0.065)] + \
#        [(0.0, 100.0), (18.6, 18.67), (1e-2, 100.0), (1e-2, 100.0), (25.85, 25.95), (8.1, 8.15), (88.85, 88.95), (99.7, 99.8), (30.6, 30.7)]


# submodel 1
# bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(1e-2, 100.0)]*len(network.get_c_graph().get_species())

print(GA.get_input_vector())

# sys.exit()

# print(network.get_c_graph().get_reactions() + network.get_c_graph().get_species())
#
# print(GA.get_decision_vector())
#

# sympy.pprint(GA.get_independent_odes_subs().jacobian(GA.get_independent_species()))
#
# sys.exit()
# params_for_global_min, obj_fun_vals = GA.run_optimization(bounds=bnds, iterations=5, seed=0, print_flag=True,
#                                                           dual_annealing_iters=1000, confidence_level_flag=True)


cons = [{'type': 'ineq', 'fun': lambda x:  x[9] - 2.0*x[8]}, {'type': 'eq', 'fun': lambda x:  x[16]}]

params_for_global_min, obj_fun_vals = GA.run_optimization(bounds=bnds, iterations=5, seed=0, print_flag=True,
                                                          dual_annealing_iters=1000, confidence_level_flag=True,
                                                          constraints=cons)


# params_for_global_min, obj_fun_vals = GA.run_optimization(bounds=bnds, iterations=5, seed=0, print_flag=True,
#                                                           confidence_level_flag=True)

# params_for_global_min, obj_fun_vals = GA.run_optimization(bounds=bnds, iterations=10, seed=0, print_flag=True,
#                                                           main_iters=1000, confidence_level_flag=True)

# sys.exit()

# params_for_global_min, obj_fun_vals, my_rank = GA.run_mpi_optimization(bounds=bnds, iterations=100, seed=0, print_flag=True,
#                                                                        dual_annealing_iters=1000, confidence_level_flag=True)

GA.generate_report()

# if my_rank == 0:
    # numpy.save('./num_cont_nuts_model_3/nuts_params.npy', params_for_global_min)
    # numpy.save('./num_cont_nuts_model/nuts_params.npy', params_for_global_min)

numpy.save('./num_cont_nuts_model/params.npy', params_for_global_min)

# params_for_global_min = numpy.load('./num_cont_nuts_model/params.npy')


# params_for_global_min = numpy.load('./num_cont_nuts_model_3/params_non_overlap.npy')

# params_for_global_min = numpy.load('./num_cont_nuts_model_3/params_non_overlap_1000.npy')


multistable_param_ind, plot_specifications = GA.run_greedy_continuity_analysis(species=response, parameters=params_for_global_min, print_lbls_flag=True,
                                                                               auto_parameters={'PrincipalContinuationParameter': signal},
                                                                               dir_path='./num_cont_nuts_model')

# multistable_param_ind, plot_specifications = GA.run_continuity_analysis(species=response, parameters=params_for_global_min, print_lbls_flag=True,
#                                                                                auto_parameters={'PrincipalContinuationParameter': signal},
#                                                                                dir_path='./num_cont_nuts_model')

sys.exit()
# cont_return_vals = [multistable_param_ind, plot_specifications]
#
# with open("./num_cont_nuts_model_3/cont_vals.dill", 'wb') as f:
#     dill.dump(cont_return_vals, f)
#
# sys.exit()

old_params = params_for_global_min[0]

# print(old_params)

lambda_vars = GA.get_variables_for_lambda_functions()

temp_array = numpy.zeros((1,len(lambda_vars)), dtype=numpy.float64)

reactions = network.get_c_graph().get_reactions()

species = network.get_c_graph().get_species()

# print(temp_array)

temp_array[0, 0:len(reactions)-2] = old_params[0:len(reactions)-2]

# filling in reaction re16
temp_array[0, len(reactions)-2] = 1.0

# filling in reaction re18
temp_array[0, len(reactions)-1] = 2.0

# print(temp_array)

indp_species = GA.get_independent_species()

# filling in old species concentrations
temp_array[0, len(reactions):len(reactions) + len(indp_species) -1] = old_params[len(reactions)-2: len(reactions)-2 + len(indp_species) -1]

# print(temp_array)

temp_array[0, len(reactions)-2 + len(indp_species)+1] = 10.0 # 0.0

# print("hi")
# print(temp_array)

# filling in conservation laws
temp_array[0, len(reactions)-2 + len(indp_species)+2] = old_params[len(reactions)-2 + len(indp_species)-1]

temp_array[0, len(reactions)-2 + len(indp_species)+3] = old_params[len(reactions)-2 + len(indp_species)]

# print(temp_array)
#
# print(lambda_vars)

params_for_global_min = [temp_array]

sympy.pprint(GA.get_independent_odes_subs())

sympy.pprint(GA.get_independent_odes_subs().jacobian(GA.get_independent_species()))


multistable_param_ind, plot_specifications = GA.run_greedy_continuity_analysis(species=response, parameters=params_for_global_min[0], print_lbls_flag=True,
                                                                               auto_parameters={'PrincipalContinuationParameter': signal},
                                                                               dir_path='./num_cont_nuts_model_3_4c')
sys.exit()

with open("./num_cont_nuts_model_3/cont_vals.dill", 'rb') as f:
    out = dill.load(f)

print(params_for_global_min)
# sys.exit()

multistable_param_ind = [0] #out[0]
print(multistable_param_ind)
plot_specifications = out[1]

print(plot_specifications)
# sys.exit()

# Parameters that produced bistability.
# re* are kinetic constants. Units can be found here help(network.get_physiological_range).
df = pandas.DataFrame(numpy.vstack([params_for_global_min[i] for i in multistable_param_ind]).T,
                      columns=["set" + str(i + 1) for i in multistable_param_ind],
                      index=[str(i) for i in GA.get_variables_for_lambda_functions()])

print(df)

################## selected parameter set #########################
decision_vector_values = numpy.array(df['set1'])
plot_specifications = plot_specifications[0]  # warning, overwriting variable!!!


################ ODEs ###################################
# print("Original ODEs")
odes = network.get_c_graph().get_ode_system()
# sympy.pprint(odes)

# why we need this? String -> Sympy objects
# construct sympy form of reactions and species
sympy_reactions = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_reactions()]
sympy_species = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_species()]
# joining together
lambda_inputs = sympy_reactions + sympy_species
# creating a lambda function for each ODE to
ode_lambda_functions = [sympy.utilities.lambdify(lambda_inputs, odes[i]) for i in range(len(odes))]

############################### kinetic constants ########################################################
kinetic_constants = numpy.array([decision_vector_values[i] for i in range(len(sympy_reactions))])

print(kinetic_constants)

################################# Computing material conservation values ############################
# combine equilibrium species' concentrations according to conservation relationships
conservation_values = numpy.array(decision_vector_values[len(sympy_reactions) + len(GA.get_independent_species()):])

print(conservation_values)

print(sympy_species)

print(GA.get_conservation_laws())



################################# starting concentrations ############################################
# this assumes that a chemical moiety in one state (species) and other species containing this moiety are zero
# assignment of conservation values to species requires exploring the model in CellDesigner
# y_fwd = [0.0, 0.0, conservation_values[0], 0.0, conservation_values[1], 0.0, 0.0]
# y_rev = [conservation_values[1], 0.0, conservation_values[0], 0.0, 0.0, 0.0, 0.0]

cons_c1_half = conservation_values[0]/2.0

y_fwd = [0.0, 0.0, cons_c1_half, 0.0, conservation_values[1], 0.0, cons_c1_half]
y_rev = [conservation_values[1], 0.0, cons_c1_half, 0.0, 0.0, 0.0, cons_c1_half]

############ simulation ###################
# computing dy/dt increments
def f(cs, t, ks, ode_lambda_func):
    return [i(*tuple(ks), *tuple(cs)) for i in ode_lambda_func]  # dy/dt

# def sim_fun_fwd(x):
#     y_fwd[2] = x  # updating s1 concentration or C3
#     return itg.odeint(f, y_fwd, t, args=(kinetic_constants, ode_lambda_functions))
#
# def sim_fun_rev(x):
#     y_rev[2] = x # updating s2 concentration
#     return itg.odeint(f, y_rev, t, args=(kinetic_constants, ode_lambda_functions))

def sim_fun_fwd(x):
    cons_c1_half = x / 2.0
    y_fwd[2] = cons_c1_half
    y_fwd[6] = cons_c1_half
    return itg.odeint(f, y_fwd, t, args=(kinetic_constants, ode_lambda_functions))

def sim_fun_rev(x):
    cons_c1_half = x / 2.0
    y_rev[2] = cons_c1_half
    y_rev[6] = cons_c1_half
    return itg.odeint(f, y_rev, t, args=(kinetic_constants, ode_lambda_functions))

# starting and ending time in seconds, number of data points
t = numpy.linspace(0.0, 60000.0, 10000)
# signal parameter scanning range and data points. Forward scan.
C3_scan = numpy.linspace(*plot_specifications[0], 60)
sim_res_fwd = [sim_fun_fwd(i) for i in C3_scan]  # occupies sys.getsizeof(sim_res_rev[0])*len(sim_res_rev)/2**20 Mb
# Reverse C3_scan. Reverse means that s2 is already high and signal is decreasing.
sim_res_rev = [sim_fun_rev(i) for i in numpy.flip(C3_scan)]

out = pandas.DataFrame(columns=['dir', 'signal', 'time'] + network.get_c_graph().get_species())
for i in range(len(sim_res_fwd)):
    out_i = pandas.DataFrame(sim_res_fwd[i], columns=out.columns[3:])
    out_i['time'] = t
    out_i['signal'] = C3_scan[i]
    out_i['dir'] = 'Low $[S^{**}]$'
    out = pandas.concat([out, out_i[out.columns]])
for i in range(len(sim_res_rev)):
    out_i = pandas.DataFrame(sim_res_rev[i], columns=out.columns[3:])
    out_i['time'] = t
    out_i['signal'] = numpy.flip(C3_scan)[i]
    out_i['dir'] = 'High $[S^{**}]$'
    out = pandas.concat([out, out_i[out.columns]])
out.to_csv("./num_cont_nuts_model_3_4c/sim2.txt", sep="\t", index=False)

###################### plotting ##################################
g = (ggplot(out, aes('time', response, group='signal', color='signal'))
     + geom_line(size=0.5)
     + ylim(0, 202)
     + labs(x="time", y="$[S^{**}]$")
     + scale_color_distiller(palette='RdYlBu', type="diverging", name="$B_{tot}$")
     + facet_wrap('~dir')
     + theme_bw())
g.save(filename="./num_cont_nuts_model_3_4c/sim_fwd_rev2.png", format="png", width=8, height=4, units='in', verbose=False)

eq = out[out.time == max(out.time)]

g = (ggplot(eq)
     + aes(x='signal', y=response, color='dir')
     + labs(x="$B_{tot}$", y="$[S^{**}]$", color="")
     + geom_path(size=2, alpha=0.5)
     + geom_point(color="black")
     + theme_bw()
     + geom_point(color="black")
     + annotate("point", x=plot_specifications[2][0][0], y=plot_specifications[2][0][1], colour="red", shape="*",
                size=3.5)
     + annotate("text", x=plot_specifications[2][0][0], y=plot_specifications[2][0][1],
                label=plot_specifications[2][0][2]))
g.save(filename="./num_cont_nuts_model_3_4c/sim_bif_diag2.png", format="png", width=6, height=4, units='in', verbose=False)
