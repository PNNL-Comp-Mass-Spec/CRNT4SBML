import sys
sys.path.insert(0, "..")

import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/Fig1Ci.xml")

network.basic_report()

network.print_c_graph()

opt = network.get_mass_conservation_approach()

bounds, concentration_bounds = opt.get_optimization_bounds()

# params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds,
#                                                                      concentration_bounds=concentration_bounds)

import numpy
#numpy.save('params.npy', params_for_global_min)
params_for_global_min = numpy.load('params.npy')

multistable_param_ind, plot_specifications = opt.run_greedy_continuity_analysis(species="s15",
                                                                                parameters=params_for_global_min[[3]],
                                                                                auto_parameters={'PrincipalContinuationParameter': 'C3'})

opt.generate_report()

import sympy

print("Original ODEs")
odes = network.get_c_graph().get_ode_system()
sympy.pprint(odes)

decision_vector_values = params_for_global_min[3]

species_concentrations = []
for i in opt.get_concentration_funs():
    species_concentrations.append(i(*tuple(decision_vector_values)))

print("species_concentrations")
print(species_concentrations)
print("species")
print(network.get_c_graph().get_species())

conservation_values = network.get_c_graph().get_b()*sympy.Matrix([species_concentrations]).T

print("conservation_values")
print(conservation_values)

print("conservation laws")
print(opt.get_conservation_laws())

# construct sympy form of reactions and species
sympy_reactions = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_reactions()]
sympy_species = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_species()]

# creating the input for each ode lambda function
lambda_inputs = sympy_reactions + sympy_species

# creating a lambda function for each ODE to
ode_lambda_functions = []
for i in range(len(odes)):
    ode_lambda_functions += [sympy.utilities.lambdify(lambda_inputs, odes[i])]

# filling in the values for the reaction rates using optimization values
input_vals = numpy.zeros(len(lambda_inputs))
for i in range(len(sympy_reactions)):
    input_vals[i] = decision_vector_values[i]

# solve the system dy/dt = f(y, t)
def f(y, t, inputs, ode_lambda_func, start_ind):

    # setting species concentrations
    inputs[start_ind:] = y

    # the model equations
    ode_vals = []
    for i in ode_lambda_func:
        ode_vals.append(i(*tuple(inputs)))

    return ode_vals


# first index of a species in input_vals
start_index = len(sympy_reactions)

# solve the ODEs
import scipy.integrate

from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
else:
    import matplotlib.pyplot as plt

fig = plt.subplots(1, 2, figsize=(15, 8))

# time grid
# final_time = 100000.0  # index 2
final_time = 2000.0 # index 3

t = numpy.linspace(0.0, final_time, 10000)

# s15_init = 2.0e5 # index 2
s15_init = 4.0e5 # index 3

# C3 values
# C3_vec = numpy.linspace(8.779e5, 8.787e5, 30) # index 2
C3_vec = numpy.linspace(8.35e5, 8.42e5, 30) # index 3

C3_min = min(C3_vec)
C3_max = max(C3_vec)

plt.subplot(1, 2, 1)

# Setting up a colormap that's a simple transtion
mymap = matplotlib.colors.LinearSegmentedColormap.from_list('mycolors', ['blue', 'red'])

CS3 = plt.cm.ScalarMappable(cmap=mymap, norm=plt.Normalize(vmin=C3_min, vmax=C3_max))

C1 = conservation_values[0]
C2 = conservation_values[1]

steady_state_vals1 = numpy.zeros(len(C3_vec))
for i in range(len(C3_vec)):
    y0 = [C3_vec[i], C2, 0.0, 0.0, C1, 0.0, 0.0]
    soln1 = scipy.integrate.odeint(f, y0, t, args=(input_vals, ode_lambda_functions, start_index))
    s15_sol = soln1[:, 6]
    steady_state_vals1[i] = soln1[-1, 6]
    r = (C3_vec[i] - C3_min) / (C3_max - C3_min)
    g = 0
    b = 1 - r
    plt.plot(t, s15_sol, color=(r, g, b))

plt.xlabel("time")
plt.ylabel("Concentration of s15")
plt.ylim(plot_specifications[0][1])
plt.ticklabel_format(axis='both', style='sci', scilimits=(-2, 2))
plt.title("Starting at s15 = 0.0")

plt.subplot(1, 2, 2)

steady_state_vals2 = numpy.zeros(len(C3_vec))
for i in range(len(C3_vec)):
    y0 = [C3_vec[i]-2*s15_init, C2, 0.0, 0.0, C1, 0.0, s15_init]
    soln1 = scipy.integrate.odeint(f, y0, t, args=(input_vals, ode_lambda_functions, start_index))
    s15_sol = soln1[:, 6]
    steady_state_vals2[i] = soln1[-1, 6]
    r = (C3_vec[i] - C3_min) / (C3_max - C3_min)
    g = 0
    b = 1 - r
    plt.plot(t, s15_sol, color=(r, g, b))

plt.xlabel("time")
plt.ylabel("Concentration of s15")
plt.ylim(plot_specifications[0][1])
plt.ticklabel_format(axis='both', style='sci', scilimits=(-2, 2))
plt.title("Starting at s15 = " + str(s15_init))

clb = plt.colorbar(CS3)
clb.set_label("C3")

plt.savefig('./num_cont_graphs/simulating_s15.png')

fig = plt.figure()

plt.plot(C3_vec, steady_state_vals1, 'bo')
plt.plot(C3_vec, steady_state_vals2, 'bo')
plt.xlabel("C3")
plt.ylabel("s15 at time " + str(final_time))
plt.ticklabel_format(axis='both', style='sci', scilimits=(-2, 2))
plt.title("Taking a Slice")

plt.savefig('./num_cont_graphs/bistability_ODEs.png')
plt.close(fig)
