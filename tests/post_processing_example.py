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

multistable_param_ind = opt.run_greedy_continuity_analysis(species="s15", parameters=params_for_global_min[[2]],
                                                           auto_parameters={'PrincipalContinuationParameter': 'C3'})

opt.generate_report()

import sympy

print("Original ODEs")
odes = network.get_c_graph().get_ode_system()
sympy.pprint(odes)

decision_vector_values = params_for_global_min[2]

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

print("decision vector")
print(decision_vector_values)

# construct sympy form of reactions and species
sympy_reactions = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_reactions()]
sympy_species = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_species()]

# creating the input for each ode lambda function
lambda_inputs = sympy_reactions + sympy_species

print("lambda_inputs")
print(lambda_inputs)

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
    ind = start_ind
    for i in y:
        inputs[ind] = i
        ind += 1

    # the model equations
    ode_vals = []
    for i in ode_lambda_func:
        ode_vals.append(i(*tuple(inputs)))

    return ode_vals


# first index of a species in
start_index = len(sympy_reactions)

# solve the ODEs
import scipy.integrate

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

# time grid
t = numpy.linspace(0.0, 10.0, 10000)

y0 = [10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0]
soln1 = scipy.integrate.odeint(f, y0, t, args=(input_vals, ode_lambda_functions, start_index))

s1_sol = soln1[:, 0]
s2_sol = soln1[:, 1]
s3_sol = soln1[:, 2]
s6_sol = soln1[:, 3]
s7_sol = soln1[:, 4]
s16_sol = soln1[:, 5]
s15_sol = soln1[:, 6]

plt.plot(t, s1_sol, label='Species s1')
plt.plot(t, s2_sol, label='Species s2')
plt.plot(t, s3_sol, label='Species s3')
plt.plot(t, s6_sol, label='Species s6')
plt.plot(t, s7_sol, label='Species s7')
plt.plot(t, s16_sol, label='Species s16')
plt.plot(t, s15_sol, label='Species s15')
plt.legend(loc='upper right')
plt.xlabel("time")
plt.ylabel("Species' Concentrations")
plt.title("Toy Values")

plt.savefig('./num_cont_graphs/toy_values.png')
plt.clf()

# time grid
t = numpy.linspace(0.0, 20000.0, 1000000)

y0 = [877900.0, 760.070288586031, 0.0, 0.0, 481700.406580797, 0.0, 0.0]
soln1 = scipy.integrate.odeint(f, y0, t, args=(input_vals, ode_lambda_functions, start_index))
s15_sol1 = soln1[:, 6]

y0 = [878400.0, 760.070288586031, 0.0, 0.0, 481700.406580797, 0.0, 0.0]
soln2 = scipy.integrate.odeint(f, y0, t, args=(input_vals, ode_lambda_functions, start_index))
s15_sol2 = soln2[:, 6]


plt.plot(t, s15_sol1, label='s1 = 877900.0')
plt.plot(t, s15_sol2, label='s1 = 878400.0')
plt.xlabel("time")
plt.ylabel("Concentration of s15")
plt.legend(loc="upper right")
plt.title("Initial Values from Optimization")

plt.savefig('./num_cont_graphs/optimization_values.png')
plt.clf()