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

# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/simple_biterminal.xml")
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/simple_biterminal_v2.xml")
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts.xml")  # No, but zero value found
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts_submodel_1.xml")  # Yes
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts_submodel_2.xml")  # No, bifurcation and limit points found and zero value found
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts_submodel_3.xml")

# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts_submodel_4c.xml")
network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts_submodel_1.xml")
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts_submodel_4d.xml")

#network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/subspace_strange.xml")

signal = "C1"
response = "s11"


# network = crnt4sbml.CRNT("../sbml_files/two_dim_tk.xml")
# signal = "C1"
# response = "s1"

network.basic_report()
network.print_c_graph()

GA = network.get_general_approach(signal=signal, response=response, fix_reactions=True)

print(GA.get_conservation_laws())
# sys.exit()
# print(GA.get_fixed_reactions())
# print(GA.get_solutions_to_fixed_reactions())
#
# sympy.pprint(GA.get_independent_odes_subs())

bnds = [(0.0, 100.0)]*len(network.get_c_graph().get_reactions()) + [(0.0, 100.0)]*len(network.get_c_graph().get_species())


# bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(1e-2, 100.0)]*(len(network.get_c_graph().get_species())-1) + [(0.0, 1e-5)]

# bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(1e-2, 100.0)]*(len(network.get_c_graph().get_species())-1) + [0.0]

# bnds = [(1,2), (3,4), (5,6), (7,8), (9,10), (11, 12), (13,14), (15,16), (17,18), (19, 20), (21,22), (23,24), (25,26), (27,28), (29,30), (31,32), 0.0]

# bnds = [(1e-2, 100.0)]*len(network.get_c_graph().get_reactions()) + [(1e-2, 100.0)]*(len(network.get_c_graph().get_species())-1) + [(0.0, 1e-12)]


print(network.get_c_graph().get_reactions() + network.get_c_graph().get_species())


# print(GA.get_input_vector())

#  [re5f, re5d, re5c, re6, re7f, re7d, re7c, re8, re16, re18, s9, s10, s2, s2s9, s11, s2s10, s1] # nuts submodel 4c
# [re5f, re5d, re5c, re6, re7f, re7d, re7c, re8,             s9, s10, s2, s2s9, s11, s2s10] # nuts submodel 1


# nuts submodel 1
# bnds = [(16.5, 17.5), (92.0, 92.5), (0.01, 0.025), (0.2, 0.25), (0.78, 0.79), (3.6, 3.7), (0.15, 0.25), (0.06, 0.065)] + \
#        [(25.85, 25.95), (8.1, 8.15), (18.6, 18.67), (88.85, 88.95), (99.7, 99.8), (30.6, 30.7)]

# bnds = [(16.5, 17.5), (92.0, 92.5), (0.01, 0.025), (0.2, 0.25), (0.78, 0.79), (3.6, 3.7), (0.15, 0.25), (0.06, 0.065), (0.0, 100.0), (0.0, 100.0)] + \
#        [(25.85, 25.95), (8.1, 8.15), (18.6, 18.67), (88.85, 88.95), (99.7, 99.8), (30.6, 30.7), (0.0, 1.0)]

# bnds = [(16.5, 17.5), (92.0, 92.5), (0.01, 0.025), (0.2, 0.25), (0.78, 0.79), (3.6, 3.7), (0.15, 0.25), (0.06, 0.065), (0.0, 100.0), (0.0, 100.0)] + \
#        [(25.85, 25.95), (8.1, 8.15), (18.6, 18.67), (88.85, 88.95), (99.7, 99.8), (30.6, 30.7), (0.0, 100.0)]

# bnds = [(10.0, 20.0), (90.0, 100.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (2.0, 4.0), (0.0, 1.0), (0.0, 1.0), (0.0, 100.0), (0.0, 100.0)] + \
#        [(20.0, 30.0), (5.0, 10.0), (15.0, 20.0), (80.0, 90.0), (90.0, 100.0), (30.0, 31.0), (0.0, 100.0)]


# bnds = [(1.6e+01, 1.7e+01), (9.15e+01, 9.3e+01), (2.0e-02, 2.1e-02), (2.1e-01, 2.3e-01), (7.8e-01, 7.9e-01), (3.6e+00, 3.7e+00),
#         (1.9e-01, 1.99e-01), (6.0e-02, 6.1e-02), (2.9e+01, 2.91e+01), (2.9e+01, 2.91e+01), (2.5e+01, 2.6e+01),
#         (8.1e+00, 8.11e+00), (1.8e+01, 1.9e+01), (8.88e+01, 8.9e+01),
#         (9.95e+01, 9.99e+01), (3.0e+01, 3.1e+01), (9.95e-05, 9.99e-05)]

# print(GA.get_input_vector())
# sys.exit()

# bnds = [16.976763, 92.360763, 0.019929, 0.218282, 0.786835, 3.685272, 0.203731, 0.062582]

# print(network.get_c_graph().get_species())
# sys.exit()

# [re5f, re5d, re5c, re6, re7f, re7d, re7c, re8]

# s10      8.116689
# s2s9    88.901835
# s11     99.785741
# s2s10   30.651950

# s2 = 138.219047 - 30.651950 - 88.901835
# s9 = 253.374303 - 8.116689 - 99.785741 - 30.651950 - 88.901835


print(GA.get_decision_vector())

params_for_global_min, obj_fun_vals = GA.run_optimization(bounds=bnds, iterations=5, seed=0, print_flag=True,
                                                          dual_annealing_iters=1000, confidence_level_flag=True)


# params_for_global_min, obj_fun_vals, my_rank = GA.run_mpi_optimization(bounds=bnds, iterations=20, seed=0, print_flag=True,
#                                                                        dual_annealing_iters=1000, confidence_level_flag=True)
#
# GA.generate_report()
#
# if my_rank == 0:
numpy.save('./num_cont_nuts_model_3/params3.npy', params_for_global_min)

# params_for_global_min = numpy.load('./num_cont_nuts_model_3/params2.npy')

# params_for_global_min = numpy.load('./num_cont_nuts_model_3/params2.npy')

# print(params_for_global_min)
#
multistable_param_ind, plot_specifications = GA.run_greedy_continuity_analysis(species=response, parameters=params_for_global_min, print_lbls_flag=True,
                                                                               auto_parameters={'PrincipalContinuationParameter': signal},
                                                                               dir_path='./num_cont_nuts_model_3')

sys.exit()

# sympy.pprint(GA.get_independent_odes_subs())
# print(params_for_global_min[5])
#
# print(GA.get_variables_for_lambda_functions())


sys.exit()
multistable_param_ind, plot_specifications = GA.run_continuity_analysis(species=response, parameters=[params_for_global_min[0]], print_lbls_flag=True,
                                                                               auto_parameters={'PrincipalContinuationParameter': signal,
                                                                                                'RL0':0.0, 'RL1':100.0, 'NMX': 1000000,
                                                                                                'ITMX': 100, 'DSMAX': 1e-2,
                                                                                                'A1': 1e10, 'ITNW': 10, 'NTST': 10, 'NCOL': 10},
                                                                               dir_path='./num_cont_nuts_model_3')

sys.exit()

# if my_rank == 0:
#     print(params_for_global_min)
# print(params_for_global_min)
# print(obj_fun_vals)

#
# GA.get_full_set_of_values(params_for_global_min)

# [re5f, re5c, re7f, re7c] [re5d, re6, re7d, re8, s9, s10, s2, s2s9, s11, s2s10]
# sys.exit()

# params_for_global_min = numpy.load('./num_cont_nuts_model_3/params.npy')
# if my_rank == 0:
#     numpy.save('./num_cont_nuts_model_3/params.npy', params_for_global_min)

print("original")
sympy.pprint(network.get_c_graph().get_ode_system())

print("")

sympy.pprint(GA.get_independent_odes_subs())

jac = GA.get_independent_odes_subs().jacobian(sympy.Matrix(GA.get_independent_species()))

print("jac of submodel 4c")
sympy.pprint(jac)

print("")

#sympy.pprint(jac.nullspace())



# det_jac = jac.det(method='lu')
# print(det_jac)
# print("")
#
# sympy_species = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_species()]
#
# # print(sympy_species)
# #
# #
# # sys.exit()
#
# zero_det_jac = det_jac.subs(GA.get_independent_species()[4], sympy.S.Zero)
#
# zero_det_jac = zero_det_jac.subs(GA.get_independent_species()[4], sympy.S.Zero)
#
# print(zero_det_jac)
# print("")


#
#
# sympy_reactions = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_reactions()]
# print(sympy_reactions)
# print("")
# zero_det_jac = det_jac.subs(sympy_reactions[8], sympy_reactions[9])
# print(zero_det_jac)

# sys.exit()
#


network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts_submodel_1.xml")

signal = "C1"
response = "s11"
# network.basic_report()
# network.print_c_graph()

GA = network.get_general_approach(signal=signal, response=response, fix_reactions=False)
jac = GA.get_independent_odes_subs().jacobian(sympy.Matrix(GA.get_independent_species()))

print("nuts submodel 1")
sympy.pprint(jac)

sys.exit()

print("")

det_jac = jac.det(method='lu')

print("nuts submodel 1")
print(det_jac)











sys.exit()

print(GA.get_independent_species())
#

# sympy.pprint(GA.get_independent_odes())
#
# print(GA.get_conservation_laws())
#
# sympy_species = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_species()]
#
# snd_eq = sympy_species[0] + sympy_species[1] - sympy.Symbol('C1', positive=True)
#
# print(snd_eq)
#
# diff_system = sympy.Matrix([[GA.get_independent_odes()[0]],[snd_eq]])
# print("")
# print("Diff system")
# sympy.pprint(diff_system)
#
# jac = diff_system.jacobian(sympy.Matrix(sympy_species))
#
# sympy.pprint(jac)
#
# print(jac.det(method='lu'))
#
# sys.exit()

multistable_param_ind, plot_specifications = GA.run_greedy_continuity_analysis(species=response, parameters=params_for_global_min, print_lbls_flag=True,
                                                                               auto_parameters={'PrincipalContinuationParameter': signal},
                                                                               dir_path='./num_cont_nuts_model_3')


# for i in multistable_param_ind:
#     print(params_for_global_min[i])
#     print(obj_fun_vals[i])

# multistable_param_ind, plot_specifications = GA.run_continuity_analysis(species=response, parameters=params_for_global_min, print_lbls_flag=True,
#                                                                                auto_parameters={'PrincipalContinuationParameter': signal},
#                                                                                dir_path='./num_cont_nuts_model_3')

# multistable_param_ind, plot_specifications = GA.run_continuity_analysis(species=response, parameters=[params_for_global_min[0]], print_lbls_flag=True,
#                                                                         auto_parameters={'PrincipalContinuationParameter': signal,
#                                                                                          'RL0': 100.0, 'RL1': 300.0, 'NMX': 1000000,
#                                                                                          'ITMX': 100, 'DSMAX': 100,
#                                                                                          'A1': 1e10, 'ITNW': 100, 'NTST': 100, 'NCOL': 100},
#                                                                         dir_path='./num_cont_nuts_model')

#multistable_param_ind

# print(GA.get_input_vector())

GA.generate_report()
sys.exit()

# numpy.save('./num_cont_nuts_model/params.npy', params_for_global_min)
# numpy.save('./num_cont_graphs/params.npy', params_for_global_min)
# params_for_global_min = numpy.load('./num_cont_nuts_model/params.npy')

# cont_return_vals = [multistable_param_ind, plot_specifications]
#
# with open("./num_cont_nuts_model/cont_vals.dill", 'wb') as f:
#     dill.dump(cont_return_vals, f)

# with open("./num_cont_nuts_sub_1_model/cont_vals.dill", 'rb') as f:
#     out = dill.load(f)

print(params_for_global_min)
# sys.exit()

multistable_param_ind = [0] #out[0]
plot_specifications = [[[138.017, 138.32], [85.92064, 127.13776], [[201.095103, 99.761304, 'LP'], [138.118, 114.452, 'LP']]]]
# [[[138.017, 138.32], [85.92064, 127.13776], [[138.219, 99.7857, 'LP'], [138.118, 114.452, 'LP']]]] # out[1]

print(plot_specifications)
# sys.exit()

# Parameters that produced bistability.
# re* are kinetic constants. Units can be found here help(network.get_physiological_range).
df = pandas.DataFrame(numpy.vstack([params_for_global_min[i] for i in multistable_param_ind]).T,
                      columns=["set" + str(i + 1) for i in multistable_param_ind],
                      index=[str(i) for i in GA.get_variables_for_lambda_functions()])

print(df)

odes = network.get_c_graph().get_ode_system()
sympy.pprint(odes)

################## selected parameter set #########################
decision_vector_values = numpy.array(df['set1'])
plot_specifications = plot_specifications[0]  # warning, overwriting variable!!!


################ ODEs ###################################
print("Original ODEs")
odes = network.get_c_graph().get_ode_system()
sympy.pprint(odes)

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

################################# starting concentrations ############################################
# this assumes that a chemical moiety in one state (species) and other species containing this moiety are zero
# assignment of conservation values to species requires exploring the model in CellDesigner
y_fwd = [0.0, 0.0, conservation_values[0], 0.0, conservation_values[1], 0.0, 0.0]
y_rev = [conservation_values[1], 0.0, conservation_values[0], 0.0, 0.0, 0.0, 0.0]

# 62.898048

# 201.0951028941771 - 62.898048

# Note, the continuation parameter C3 (first position) will be varied during simulations

############ simulation ###################
# computing dy/dt increments
def f(cs, t, ks, ode_lambda_func):
    return [i(*tuple(ks), *tuple(cs)) for i in ode_lambda_func]  # dy/dt

def sim_fun_fwd(x):
    y_fwd[2] = x  # updating s1 concentration or C3
    return itg.odeint(f, y_fwd, t, args=(kinetic_constants, ode_lambda_functions))

def sim_fun_rev(x):
    y_rev[2] = x # updating s2 concentration
    return itg.odeint(f, y_rev, t, args=(kinetic_constants, ode_lambda_functions))

# starting and ending time in seconds, number of data points
t = numpy.linspace(0.0, 7000.0, 100)
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
out.to_csv("./num_cont_graphs/sim2.txt", sep="\t", index=False)

###################### plotting ##################################
g = (ggplot(out, aes('time', response, group='signal', color='signal'))
     + geom_line(size=0.5)
     + ylim(0, 202)
     + labs(x="time", y="$[S^{**}]$")
     + scale_color_distiller(palette='RdYlBu', type="diverging", name="$B_{tot}$")
     + facet_wrap('~dir')
     + theme_bw())
g.save(filename="./num_cont_graphs/sim_fwd_rev2.png", format="png", width=8, height=4, units='in', verbose=False)

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
     # + annotate("point", x=plot_specifications[2][1][0], y=plot_specifications[2][1][1], colour="red", shape="*",
     #            size=3.5)
     # + annotate("text", x=plot_specifications[2][1][0], y=plot_specifications[2][1][1],
     #            label=plot_specifications[2][1][2]))
g.save(filename="./num_cont_graphs/sim_bif_diag2.png", format="png", width=6, height=4, units='in', verbose=False)
