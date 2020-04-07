import sys
sys.path.insert(0, "..")
import crnt4sbml
import numpy
import sympy

network = crnt4sbml.CRNT("../sbml_files/insulin_signaling_motifs/Nuts_submodel_1.xml")

signal = "C1"
response = "s11"

# network.basic_report()

GA = network.get_general_approach() 

# print(GA.get_conservation_laws())
# print(network.get_c_graph().get_species())

GA.initialize_general_approach(signal=signal, response=response)

sympy_reactions = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_reactions()]
sympy_species = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_species()]
# joining together
lambda_inputs = sympy_reactions + sympy_species
# creating a lambda function for each ODE to
odes = network.get_c_graph().get_ode_system()
ode_lambda_functions = [sympy.utilities.lambdify(lambda_inputs, odes[i]) for i in range(len(odes))]
# sympy.pprint(odes)
odes_jac = odes.jacobian(sympy_species)
# sympy.pprint(odes_jac)

odes_jac_lambda = sympy.utilities.lambdify(lambda_inputs, odes_jac)

params_for_global_min = numpy.load('./num_cont_direct/params.npy')

result_x = params_for_global_min[0]


__R = len(network.get_c_graph().get_reactions())

__N = len(network.get_c_graph().get_species())

# jac_out = odes_jac_lambda(*tuple(result_x[0:__R + __N]))
#
# print(jac_out.shape[0])
#
# list_out = [list(jac_out[i]) for i in range(jac_out.shape[0])]
#
# print(list_out)


# print(params_for_global_min)

__cons_laws_sympy_lamb = GA.get_lambda_conservation_laws()

conservation_vals = [__cons_laws_sympy_lamb[i](*tuple(result_x[__R:__R + __N]))
                                 for i in range(len(__cons_laws_sympy_lamb))]

# print(conservation_vals)


def ff(t, cs, ks, ode_lambda_func):

    return [i(*tuple(ks), *tuple(cs)) for i in ode_lambda_func]

def jac_ff(t, cs, ks, odes_jac_lam):

    jac_out = odes_jac_lam(*tuple(ks), *tuple(cs))

    return [list(jac_out[i]) for i in range(jac_out.shape[0])]

initial_species_values = [0.0 for i in range(__N)]

initial_species_values[2] = conservation_vals[0]

initial_species_values[4] = conservation_vals[1]

print(f"intial species = {initial_species_values}")

# from scipy.optimize import fsolve
# eq = fsolve(f, initial_species_values, args=(0, result_x[0:__R], ode_lambda_functions))
# print(eq)

import scipy.integrate as itg
import pandas
from plotnine import ggplot, aes, geom_line, ylim, scale_color_distiller, facet_wrap, theme_bw, geom_path, geom_point, labs, annotate
from matplotlib import rc
rc('text', usetex=True)


def steady_state_finder(initial_species_values, ff, result_x, ode_lambda_functions, spec_index):


    out = itg.solve_ivp(ff, [0.0, 100.0], initial_species_values, args=(result_x[0:__R], ode_lambda_functions), method='LSODA')
    y0 = out.y[:, -1]

    i = 1
    while abs(out.y[spec_index, -1] - out.y[spec_index, 0]) > 1e-10 and i < 1000:
        tspan = [0.0 + i*100.0, 100.0 + i*100.0]
        out = itg.solve_ivp(ff, tspan, y0, args=(result_x[0:__R], ode_lambda_functions), method='LSODA')
        y0 = out.y[:, -1]
        i += 1
    if i >= 1000:
        print("Iterations greater than 1000 occured in steady state finder.")
    return out.y[:, -1]

spec_index = network.get_c_graph().get_species().index(response)    # 4
spec_index_sim = spec_index

#initial_species_values[2] = conservation_vals[0]

initial_species_values[4] = conservation_vals[1]

pcp_scan = numpy.linspace(90.0, 93.5, 60)

forward_scan = []
for i in pcp_scan:
    initial_species_values[2] = i
    steady_state = steady_state_finder(initial_species_values, ff, result_x, ode_lambda_functions, spec_index)
    forward_scan.append(steady_state[spec_index])

print(f"forward scan = {forward_scan}")

# sys.exit()

initial_species_values[4] = 0.0
initial_species_values[0] = conservation_vals[1]
# spec_index_sim = 0

reverse_scan = []
for i in pcp_scan:
    initial_species_values[2] = i
    steady_state = steady_state_finder(initial_species_values, ff, result_x, ode_lambda_functions, spec_index)
    reverse_scan.append(steady_state[spec_index])

print(f"reverse scan = {reverse_scan}")


def plot_direct_simulation(response_species, pcp_scan, forward_scan, reverse_scan, path):

    out = pandas.DataFrame(columns=['dir', 'signal'] + [response_species])
    for i in range(len(forward_scan)):
        out_i = pandas.DataFrame([forward_scan[i]], columns=[out.columns[2]])
        out_i['signal'] = pcp_scan[i]
        out_i['dir'] = 'Low concentration'
        out = pandas.concat([out, out_i[out.columns]])
    for i in range(len(reverse_scan)):
        out_i = pandas.DataFrame([reverse_scan[i]], columns=[out.columns[2]])
        out_i['signal'] = pcp_scan[i]
        out_i['dir'] = 'High concentration'
        out = pandas.concat([out, out_i[out.columns]])
    #out.to_csv(path + "/sim2.txt", sep="\t", index=False)

    g = (ggplot(out)
         + aes(x='signal', y=response_species, color='dir')
         + labs(x="$B_{tot}$", y="$[S^{**}]$", color="")
         + geom_path(size=2, alpha=0.5)
         + geom_point(color="black")
         + theme_bw()
         + geom_point(color="black"))
         # + annotate("point", x=plot_specifications[2][0][0], y=plot_specifications[2][0][1], colour="red", shape="*",
         #            size=3.5)
         # + annotate("text", x=plot_specifications[2][0][0], y=plot_specifications[2][0][1],
         #            label=plot_specifications[2][0][2]))
    # + annotate("point", x=plot_specifications[2][1][0], y=plot_specifications[2][1][1], colour="red", shape="*",
    #            size=3.5)
    # + annotate("text", x=plot_specifications[2][1][0], y=plot_specifications[2][1][1],
    #            label=plot_specifications[2][1][2]))
    g.save(filename=path+"/sim_bif_diag.png", format="png", width=6, height=4, units='in', verbose=False)


response_species = network.get_c_graph().get_species()[spec_index]
path = './num_cont_direct'
plot_direct_simulation(response_species, pcp_scan, forward_scan, reverse_scan, path)

sys.exit()

print(out[-1:])

t = numpy.linspace(0.0, 100.0, 100)
out = itg.odeint(f, initial_species_values, t, args=(result_x[0:__R], ode_lambda_functions))
temp_out = out[-1:]

print(temp_out)

t = numpy.linspace(100.0, 1000.0, 100)
out = itg.odeint(f, list(temp_out[0]), t, args=(result_x[0:__R], ode_lambda_functions))
print(out[-1:])


sys.exit()

from scipy.integrate import ode
import warnings

# backend = 'dopri5'

# solver = ode(ff, jac_ff).set_integrator(backend, nsteps=1)
# solver.set_initial_value(initial_species_values, 0.0).set_f_params(*(result_x[0:__R], ode_lambda_functions)).set_jac_params(*(result_x[0:__R], odes_jac_lambda))
#
# # solver = ode(ff).set_integrator(backend, nsteps=1)
# # solver.set_initial_value(initial_species_values, 0.0).set_f_params(*(result_x[0:__R], ode_lambda_functions))
#
# solver._integrator.iwork[2] = -1
#
# import time
#
# t1 = 100.0
# sol = []
# warnings.filterwarnings("ignore", category=UserWarning)
# start = time.time()
# while solver.t < t1:
#     #print(solver.t)
#     solver.integrate(t1, step=True)
# end = time.time()
#
# print(f"elapsed time {end - start}")
#
# warnings.resetwarnings()
#
# print(solver.t)
# print(solver.y)
