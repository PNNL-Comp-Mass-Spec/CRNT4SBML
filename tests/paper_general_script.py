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

network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/a_b.xml") # yes 10
signal = "C1"
#response = "s6"
response = "s5"
iters = 10

network.basic_report()
network.print_c_graph()

GA = network.get_general_approach(signal=signal, response=response)

print(GA.get_conservation_laws())
print(GA.get_fixed_reactions())
print(GA.get_solutions_to_fixed_reactions())


bnds = [(1e-2, 1e2)]*len(GA.get_input_vector())

# params_for_global_min, obj_fun_vals = GA.run_optimization(bounds=bnds, iterations=iters, seed=0, print_flag=False,
#                                            dual_annealing_iters=1000, confidence_level_flag=True)


#numpy.save('./basic_example_paper/params.npy', params_for_global_min)
params_for_global_min = numpy.load('./basic_example_paper/params.npy')

# multistable_param_ind, plot_specifications = GA.run_greedy_continuity_analysis(species=response, parameters=params_for_global_min,
#                                                                                auto_parameters={'PrincipalContinuationParameter': signal},
#                                                                                dir_path="./basic_example_paper")

# cont_return_vals = [multistable_param_ind, plot_specifications]
#
# with open("./basic_example_paper/cont_vals.dill", 'wb') as f:
#     dill.dump(cont_return_vals, f)

with open("./basic_example_paper/cont_vals.dill", 'rb') as f:
    out = dill.load(f)

multistable_param_ind = out[0]
plot_specifications = out[1]

#GA.generate_report()

# Parameters that produced bistability.
# re* are kinetic constants. Units can be found here help(network.get_physiological_range).
df = pandas.DataFrame(numpy.vstack([params_for_global_min[i] for i in multistable_param_ind]).T,
                      columns=["set" + str(i + 1) for i in multistable_param_ind],
                      index=[str(i) for i in GA.get_variables_for_lambda_functions()])

print(df)
sys.exit()

################## selected parameter set #########################
decision_vector_values = numpy.array(df['set2'])
plot_specifications = plot_specifications[1]  # warning, overwriting variable!!!

################ ODEs ###################################
print("Original ODEs")
odes = network.get_c_graph().get_ode_system()

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

################################# Computing material conservation values ############################
# combine equilibrium species' concentrations according to conservation relationships
conservation_values = numpy.array(decision_vector_values[len(sympy_reactions) + len(GA.get_independent_species()):])

################################# starting concentrations ############################################
# this assumes that a chemical moiety in one state (species) and other species containing this moiety are zero
# assignment of conservation values to species requires exploring the model in CellDesigner
y_fwd = [0.0, conservation_values[0], 0.0]  # [1.343145, conservation_values[0], 5.697031]
y_rev = [50.0, conservation_values[0], 0.0]  # [conservation_values[0], 9.962916618928801, 5.697031]

# Note, the continuation parameter C3 (first position) will be varied during simulations

############ simulation ###################
# computing dy/dt increments
def f(cs, t, ks, ode_lambda_func):
    return [i(*tuple(ks), *tuple(cs)) for i in ode_lambda_func]  # dy/dt

def sim_fun_fwd(x):
    y_fwd[1] = x  # updating s1 concentration or C3
    return itg.odeint(f, y_fwd, t, args=(kinetic_constants, ode_lambda_functions))

def sim_fun_rev(x):
    y_rev[1] = x # updating s2 concentration
    return itg.odeint(f, y_rev, t, args=(kinetic_constants, ode_lambda_functions))

# starting and ending time in seconds, number of data points
t = numpy.linspace(0.0, 20000.0, 10000)
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
    out_i['dir'] = 'Low [A]'
    out = pandas.concat([out, out_i[out.columns]])
for i in range(len(sim_res_rev)):
    out_i = pandas.DataFrame(sim_res_rev[i], columns=out.columns[3:])
    out_i['time'] = t
    out_i['signal'] = numpy.flip(C3_scan)[i]
    out_i['dir'] = 'High [A]'
    out = pandas.concat([out, out_i[out.columns]])
out.to_csv("./basic_example_paper/sim.txt", sep="\t", index=False)

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
     + geom_point(color="black")
     + annotate("point", x=plot_specifications[2][0][0], y=plot_specifications[2][0][1], colour="red", shape="*", size=3.5)
     + annotate("text", x=plot_specifications[2][0][0]+0.02, y=plot_specifications[2][0][1], label=plot_specifications[2][0][2])
     + annotate("point", x=plot_specifications[2][1][0], y=plot_specifications[2][1][1], colour="red", shape="*", size=3.5)
     + annotate("text", x=plot_specifications[2][1][0] - 0.02, y=plot_specifications[2][1][1], label=plot_specifications[2][1][2])
     )

g.save(filename="./basic_example_paper/sim_bif_diag.png", format="png", width=6, height=4, units='in', verbose=False)