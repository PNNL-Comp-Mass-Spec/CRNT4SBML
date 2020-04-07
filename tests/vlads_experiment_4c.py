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

# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts_submodel_4c.xml")
# network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts_submodel_1.xml")
network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts.xml")

signal = "C2" # "C1"
response = "s11"

# network.basic_report()
# network.print_c_graph()

GA = network.get_general_approach(signal=signal, response=response, fix_reactions=False)

# GA = network.get_general_approach()
# print(GA.get_conservation_laws())
# GA.initialize_general_approach(signal=signal, response=response)

print(GA.get_conservation_laws())

# sys.exit()

# submodel 4c bounds
# bnds = [(16.5, 17.5), (92.0, 92.5), (0.01, 0.025), (0.2, 0.25), (0.78, 0.79), (3.6, 3.7), (0.15, 0.25), (0.06, 0.065),    (0.9, 1.1), (1.9, 2.1)] + \
#        [(25.85, 25.95), (8.1, 8.15), (18.6, 18.67), (88.85, 88.95), (99.7, 99.8), (30.6, 30.7), (0.0, 1.0)]

# Nuts bounds
# bnds = [(1e-2, 100.0), (1e-2, 100.0), (1e-2, 100.0), (1e-2, 100.0), (0.9, 1.1), (1.9, 2.1), (16.5, 17.5), (92.0, 92.5), (0.01, 0.025), (0.2, 0.25), (0.78, 0.79), (3.6, 3.7), (0.15, 0.25), (0.06, 0.065)] + \
#        [(0.0, 1.0), (18.6, 18.67), (0.0, 1.0), (0.0, 1.0), (25.85, 25.95), (8.1, 8.15), (88.85, 88.95), (99.7, 99.8), (30.6, 30.7)]

# bnds = [(2.4, 2.42), (27.5, 28.1), (2.0, 2.15), (48.25, 48.4), (0.5, 1.1), (1.8, 2.1), (17.0, 17.5), (92.4, 92.6), (0.01, 0.025), (0.2, 0.25), (0.78, 0.79), (3.6, 3.7), (0.15, 0.25), (0.06, 0.065)] + \
#        [(0.0, 1e-12), (18.0, 18.5), (0.0, 1e-12), (0.0, 1e-12), (27.0, 27.1), (8.2, 8.3), (90.0, 90.1), (97.5, 97.9), (30.0, 30.1)]

bnds = [(2.4, 2.42), (27.5, 28.1), (2.0, 2.15), (48.25, 48.4), (0.5, 1.1), (1.8, 2.1), (17.0, 17.5), (92.4, 92.6), (0.01, 0.025), (0.2, 0.25), (0.78, 0.79), (3.6, 3.7), (0.15, 0.25), (0.06, 0.065)] + \
       [(1e-2, 100.0), (18.0, 18.5), (0.0, 0.45), (0.5, 1.0), (27.0, 27.1), (8.2, 8.3), (90.0, 90.1), (97.5, 97.9), (30.0, 30.1)]

# bnds = [(2.4, 2.42), (27.5, 28.1), (2.0, 2.15), (48.25, 48.4), (0.5, 1.1), (1.8, 2.1), (17.0, 17.5), (92.4, 92.6), (0.01, 0.025), (0.2, 0.25), (0.78, 0.79), (3.6, 3.7), (0.15, 0.25), (0.06, 0.065)] + \
#        [(0.0, 100.0), (18.0, 18.5), (0.5, 1.0), (0.0, 0.45), (27.0, 27.1), (8.2, 8.3), (90.0, 90.1), (97.5, 97.9), (30.0, 30.1)]

# bnds = [(2.4, 2.42), (27.5, 28.1), (2.0, 2.15), (48.25, 48.4), (0.5, 1.1), (1.8, 2.1), (17.0, 17.5), (92.4, 92.6), (0.01, 0.025), (0.2, 0.25), (0.78, 0.79), (3.6, 3.7), (0.15, 0.25), (0.06, 0.065)] + \
#        [(0.0, 100.0), (18.0, 18.5), (0.5, 1.0), (0.0, 0.45), (27.0, 27.1), (8.2, 8.3), (90.0, 90.1), (97.5, 97.9), (30.0, 30.1)]

print(GA.get_input_vector())

# re3 != re4*s6/s3
# re7c != re8*s11/s2s10

# 0.06497873603238077*71.82806697466788/31.01894105572572

# sys.exit()

# print(network.get_c_graph().get_reactions() + network.get_c_graph().get_species())
#
# print(GA.get_decision_vector())
print(GA.get_fixed_reactions())
print(GA.get_solutions_to_fixed_reactions())

# sys.exit()

#
# params_for_global_min, obj_fun_vals = GA.run_optimization(bounds=bnds, iterations=5, seed=0, print_flag=True,
#                                                           dual_annealing_iters=1000, confidence_level_flag=True)
#
# GA.generate_report()
#
# sys.exit()

params_for_global_min, obj_fun_vals, my_rank = GA.run_mpi_optimization(bounds=bnds, iterations=100, seed=0, print_flag=True,
                                                                       dual_annealing_iters=1000, confidence_level_flag=True)

GA.generate_report()

if my_rank == 0:
    # numpy.save('./num_cont_nuts_model_3/nuts_params.npy', params_for_global_min)
    # numpy.save('./num_cont_nuts_model_new/nuts_params.npy', params_for_global_min)

    # numpy.save('./num_cont_nuts_model_new/nuts_params_2.npy', params_for_global_min)
    # numpy.save('./num_cont_nuts_model_new/nuts_params_3.npy', params_for_global_min)
    numpy.save('./num_cont_nuts_model_new/nuts_params_4.npy', params_for_global_min)

# numpy.save('./num_cont_nuts_model_3/params_non_overlap.npy', params_for_global_min)


# params_for_global_min = numpy.load('./num_cont_nuts_model_3/params_non_overlap.npy')

# params_for_global_min = numpy.load('./num_cont_nuts_model_3/params_non_overlap_1000.npy')

# params_for_global_min = numpy.load('./num_cont_nuts_model_new/nuts_params.npy')

# params_for_global_min = numpy.load('./num_cont_nuts_model_new/nuts_params_2.npy')

# params_for_global_min = numpy.load('./num_cont_nuts_model_new/nuts_params_3.npy')

# params_for_global_min = numpy.load('./num_cont_nuts_model_new/nuts_params_4.npy')

# print(params_for_global_min)

sys.exit()

multistable_param_ind, plot_specifications = GA.run_greedy_continuity_analysis(species=response, parameters=params_for_global_min, print_lbls_flag=True,
                                                                               auto_parameters={'PrincipalContinuationParameter': signal},
                                                                               dir_path='./num_cont_nuts_model_new')
sys.exit()
# lambda_funcs = GA.get_indpendent_odes_lambda_function()
#
# jacobian = GA.get_independent_odes_subs().jacobian(GA.get_independent_species())
#
# sympy.pprint(jacobian)

# temp_jac = jacobian[:, :]

# for i in range(jacobian.shape[0]):
#     for j in range(jacobian.shape[1]):
#         for ii in range(len(GA.get_fixed_reactions())):
#             if ii not in [2, 6]:
#                 jacobian[i, j] = jacobian[i, j].subs(GA.get_fixed_reactions()[ii], GA.get_solutions_to_fixed_reactions()[ii])
#         jacobian[i, j] = sympy.simplify(jacobian[i, j])

# sympy.pprint(jacobian)
# print(jacobian)
# sys.exit()



# det_jac = jacobian.det(method='lu')

# print(det_jac)



# print(GA.get_determinant_of_jacobian())

# det_jac = GA.get_determinant_of_jacobian()
#
# print(det_jac)

# sys.exit()

# temp_jac = jacobian[:, :]
#
# for i in range(len(GA.get_fixed_reactions())):
#     temp_jac = temp_jac.subs(GA.get_fixed_reactions()[i], GA.get_solutions_to_fixed_reactions()[i])
#
#
# sympy.pprint(temp_jac)
#
# sys.exit()

# rref_jac = temp_jac.rref()
#
# sympy.pprint(rref_jac)

# print(f"sympy jacobian rank {jacobian.rank()}")

# lamb_jac = sympy.utilities.lambdify(GA.get_variables_for_lambda_functions(), jacobian)
#
# array_val = params_for_global_min[0][:]
#
# print(array_val)
# print(f"array val = {array_val[21]}")
#
# array_val[21] = 300.0
#
# temp = lamb_jac(*tuple(array_val))
#
# from numpy.linalg import matrix_rank
#
# print(f"matrix rank = {matrix_rank(temp)}")
#
# print(numpy.linalg.svd(temp)[1])
#
# print("sympy array")
#
# subs_vals = sympy.Matrix(temp)
#
# sympy.pprint(subs_vals)
#
# sympy.pprint(subs_vals.rref())


# for i in range(temp.shape[0]):
#     print(temp[i,:])

# print(GA.get_variables_for_lambda_functions())
#
# for i in lambda_funcs:
#     print(i(*tuple(array_val)))
#
# sys.exit()


# multistable_param_ind, plot_specifications = GA.run_greedy_continuity_analysis(species=response, parameters=params_for_global_min, print_lbls_flag=True,
#                                                                         auto_parameters={'PrincipalContinuationParameter': signal},
#                                                                         dir_path='./num_cont_nuts_model_new')
#
# sys.exit()

# cont_return_vals = [multistable_param_ind, plot_specifications]
#
# with open("./num_cont_nuts_model_3/cont_vals.dill", 'wb') as f:
#     dill.dump(cont_return_vals, f)
#
# sys.exit()

# old_params = params_for_global_min[0]

# print(old_params)

# lambda_vars = GA.get_variables_for_lambda_functions()
#
# temp_array = numpy.zeros((1,len(lambda_vars)), dtype=numpy.float64)
#
# reactions = network.get_c_graph().get_reactions()
#
# species = network.get_c_graph().get_species()

# print(temp_array)

# temp_array[0, 0:len(reactions)-2] = old_params[0:len(reactions)-2]

# filling in reaction re16
# temp_array[0, len(reactions)-2] = 1.0

# filling in reaction re18
# temp_array[0, len(reactions)-1] = 2.0

# print(temp_array)

# indp_species = GA.get_independent_species()

# filling in old species concentrations
# temp_array[0, len(reactions):len(reactions) + len(indp_species) -1] = old_params[len(reactions)-2: len(reactions)-2 + len(indp_species) -1]

# print(temp_array)

# temp_array[0, len(reactions)-2 + len(indp_species)+1] = 10.0 # 0.0

# print("hi")
# print(temp_array)

# filling in conservation laws
# temp_array[0, len(reactions)-2 + len(indp_species)+2] = old_params[len(reactions)-2 + len(indp_species)-1]
#
# temp_array[0, len(reactions)-2 + len(indp_species)+3] = old_params[len(reactions)-2 + len(indp_species)]

# print(temp_array)
#
# print(lambda_vars)

# params_for_global_min = [temp_array]
#
# sympy.pprint(GA.get_independent_odes_subs())
#
# sympy.pprint(GA.get_independent_odes_subs().jacobian(GA.get_independent_species()))
#
#
# multistable_param_ind, plot_specifications = GA.run_greedy_continuity_analysis(species=response, parameters=params_for_global_min[0], print_lbls_flag=True,
#                                                                                auto_parameters={'PrincipalContinuationParameter': signal},
#                                                                                dir_path='./num_cont_nuts_model_3_4c')
# sys.exit()

with open("./num_cont_nuts_model_3/cont_vals.dill", 'rb') as f:
    out = dill.load(f)

# print(params_for_global_min)
# sys.exit()

multistable_param_ind = [i for i in range(len(params_for_global_min))] # out[0]
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

# re2 = re4*s6/s3
# 2.068267 = 1.917226*0.440630/0.910597
# re7c = re8*s11/s2s10
# 0.207868 = 0.063727*97.885203/30.009321

# 0.060036*97.051925/30.082454

sys.exit()
# ################## selected parameter set #########################
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

sys.exit()

################################# starting concentrations ############################################
# this assumes that a chemical moiety in one state (species) and other species containing this moiety are zero
# assignment of conservation values to species requires exploring the model in CellDesigner
# y_fwd = [0.0, 0.0, conservation_values[0], 0.0, conservation_values[1], 0.0, 0.0]
# y_rev = [conservation_values[1], 0.0, conservation_values[0], 0.0, 0.0, 0.0, 0.0]

# [s9, s10, s2, s2s9, s11, s2s10]
# C1 = 1.0*s2 + 1.0*s2s10 + 1.0*s2s9
# C2 = 1.0*s10 + 1.0*s11 + 1.0*s2s10 + 1.0*s2s9 + 1.0*s9


y_fwd = [0.0, conservation_values[1], 0.0, 0.0, 0.0, 0.0, 0.0, conservation_values[0], 0.0]
y_rev = [0.0, conservation_values[1], 0.0, 0.0, conservation_values[0], 0.0, 0.0, 0.0, 0.0]

# y_fwd = [conservation_values[1], 0.0, 0.0, 0.0, conservation_values[0], 0.0, 0.0, 0.0, 0.0]
# y_rev = [conservation_values[1], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, conservation_values[0], 0.0]

# cons_c1_half = conservation_values[0]/2.0
#
# y_fwd = [0.0, 0.0, cons_c1_half, 0.0, conservation_values[1], 0.0, cons_c1_half]
# y_rev = [conservation_values[1], 0.0, cons_c1_half, 0.0, 0.0, 0.0, cons_c1_half]

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
    y_fwd[1] = x  # updating s1 concentration or C3
    return itg.odeint(f, y_fwd, t, args=(kinetic_constants, ode_lambda_functions))

def sim_fun_rev(x):
    y_rev[1] = x # updating s2 concentration
    return itg.odeint(f, y_rev, t, args=(kinetic_constants, ode_lambda_functions))

# def sim_fun_fwd(x):
#     cons_c1_half = x / 2.0
#     y_fwd[2] = cons_c1_half
#     y_fwd[6] = cons_c1_half
#     return itg.odeint(f, y_fwd, t, args=(kinetic_constants, ode_lambda_functions))
#
# def sim_fun_rev(x):
#     cons_c1_half = x / 2.0
#     y_rev[2] = cons_c1_half
#     y_rev[6] = cons_c1_half
#     return itg.odeint(f, y_rev, t, args=(kinetic_constants, ode_lambda_functions))

# starting and ending time in seconds, number of data points
t = numpy.linspace(0.0, 70000.0, 10000)
# signal parameter scanning range and data points. Forward scan.
# C3_scan = numpy.linspace(138.017, 138.32, 100)
# C3_scan = numpy.linspace(139.0, 141.0, 100)
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
out.to_csv("./num_cont_nuts_model/sim2.txt", sep="\t", index=False)

###################### plotting ##################################
g = (ggplot(out, aes('time', response, group='signal', color='signal'))
     + geom_line(size=0.5)
     # + ylim(0, 202)
     + labs(x="time", y="$[S^{**}]$")
     + scale_color_distiller(palette='RdYlBu', type="diverging", name="$B_{tot}$")
     + facet_wrap('~dir')
     + theme_bw())
g.save(filename="./num_cont_nuts_model/sim_fwd_rev2.png", format="png", width=8, height=4, units='in', verbose=False)

eq = out[out.time == max(out.time)]

# with pandas.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
#     print(eq)

g = (ggplot(eq)
     + aes(x='signal', y=response, color='dir')
     + labs(x="$B_{tot}$", y="$[S^{**}]$", color="")
     + geom_path(size=2, alpha=0.5)
     + geom_point(color="black")
     + theme_bw()
     + geom_point(color="black"))
     # + annotate("point", x=plot_specifications[2][0][0], y=plot_specifications[2][0][1], colour="red", shape="*",
     #            size=3.5)
     # + annotate("text", x=plot_specifications[2][0][0], y=plot_specifications[2][0][1],
     #            label=plot_specifications[2][0][2]))
g.save(filename="./num_cont_nuts_model/sim_bif_diag2.png", format="png", width=6, height=4, units='in', verbose=False)
