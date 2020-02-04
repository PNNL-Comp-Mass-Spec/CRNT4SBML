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

network = crnt4sbml.CRNT("../sbml_files/Song.xml") # yes for 100
signal = "C1"
response = "s2"
iters = 100

network.basic_report()
network.print_c_graph()

GA = network.get_general_approach(signal=signal, response=response)

print(GA.get_conservation_laws())
print(GA.get_fixed_reactions())
print(GA.get_solutions_to_fixed_reactions())

print(network.get_c_graph().get_reactions())
print(network.get_c_graph().get_species())
print(GA.get_decision_vector())

print(network.get_c_graph().get_s().rank())

sympy.pprint(GA.get_independent_species())

#sys.exit()
#bnds = GA.get_optimization_bounds()
bnds = [(1e-3, 6.0)]*len(network.get_c_graph().get_reactions()) + [(1e-3, 1000.0)]*len(network.get_c_graph().get_species())

# params_for_global_min, obj_fun_vals = GA.run_optimization(bounds=bnds, iterations=iters, seed=0, print_flag=False,
#                                            dual_annealing_iters=1000, confidence_level_flag=True)


# numpy.save('./num_cont_song_example/params.npy', params_for_global_min)
params_for_global_min = numpy.load('./num_cont_song_example/params.npy')

# multistable_param_ind, plot_specifications = GA.run_greedy_continuity_analysis(species=response, parameters=params_for_global_min,
#                                                                                auto_parameters={'PrincipalContinuationParameter': signal},
#                                                                                dir_path="./num_cont_song_example")

# cont_return_vals = [multistable_param_ind, plot_specifications]
#
# with open("./num_cont_song_example/cont_vals.dill", 'wb') as f:
#     dill.dump(cont_return_vals, f)

with open("./num_cont_song_example/cont_vals.dill", 'rb') as f:
    out = dill.load(f)

multistable_param_ind = out[0]
plot_specifications = out[1]

#GA.generate_report()

#sys.exit()

# Parameters that produced bistability.
# re* are kinetic constants. Units can be found here help(network.get_physiological_range).
df = pandas.DataFrame(numpy.vstack([params_for_global_min[i] for i in multistable_param_ind]).T,
                      columns=["set" + str(i + 1) for i in multistable_param_ind],
                      index=[str(i) for i in GA.get_variables_for_lambda_functions()])

print(df)

sys.exit()

################## selected parameter set #########################
decision_vector_values = numpy.array(df['set2'])
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
y_fwd = [conservation_values[1], conservation_values[0], 0.0, 0.0, 0.0, 0.0]
y_rev = [0.0, conservation_values[0], 0.0, conservation_values[1], 0.0, 0.0]

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
t = numpy.linspace(0.0, 2000.0, 3000)
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
    out_i['dir'] = 'Low $[S_p]$'
    out = pandas.concat([out, out_i[out.columns]])
for i in range(len(sim_res_rev)):
    out_i = pandas.DataFrame(sim_res_rev[i], columns=out.columns[3:])
    out_i['time'] = t
    out_i['signal'] = numpy.flip(C3_scan)[i]
    out_i['dir'] = 'High $[S_p]$'
    out = pandas.concat([out, out_i[out.columns]])
out.to_csv("./num_cont_song_example/sim.txt", sep="\t", index=False)

###################### plotting ##################################
g = (ggplot(out, aes('time', response, group='signal', color='signal'))
     + geom_line(size=0.5)
     + ylim(0, 1000)
     + labs(x="time", y="$[S_p]$")
     + scale_color_distiller(palette='RdYlBu', type="diverging", name="$K_{tot}$")
     + facet_wrap('~dir')
     + theme_bw())
g.save(filename="./num_cont_song_example/sim_fwd_rev.png", format="png", width=8, height=4, units='in', verbose=False)

eq = out[out.time == max(out.time)]
g = (ggplot(eq)
     + aes(x='signal', y=response, color='dir')
     + labs(x="$K_{tot}$", y="$[S_p]$", color="")
     + geom_path(size=2, alpha=0.5)
     + geom_point(color="black")
     + theme_bw()
     + geom_point(color="black")
     + annotate("point", x=plot_specifications[2][0][0], y=plot_specifications[2][0][1], colour="red", shape="*",
                size=3.5)
     + annotate("text", x=plot_specifications[2][0][0] + 5.00, y=plot_specifications[2][0][1],
                label=plot_specifications[2][0][2])
     + annotate("point", x=plot_specifications[2][1][0], y=plot_specifications[2][1][1], colour="red", shape="*",
                size=3.5)
     + annotate("text", x=plot_specifications[2][1][0] - 5.00, y=plot_specifications[2][1][1],
                label=plot_specifications[2][1][2]))
g.save(filename="./num_cont_song_example/sim_bif_diag.png", format="png", width=6, height=4, units='in', verbose=False)