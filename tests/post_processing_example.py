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

odes = opt.get_independent_odes()

import sympy

print("Original ODEs")
sympy.pprint(network.get_c_graph().get_ode_system())
print("")
print("Independent ODEs")
sympy.pprint(odes)

species = opt.get_independent_species()

print("Independent Species")
print(species)

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

#creating the input for each ode lambda function
lambda_inputs = [sympy.Symbol(i, positive=True) for i in network.get_c_graph().get_reactions()] + species + \
                [sympy.Symbol('C1', real=True), sympy.Symbol('C2', real=True), sympy.Symbol('C3', real=True)]

print(lambda_inputs)

#creating a lambda function for each ODE to
ode_lambda_functions = []
for i in range(len(odes)):
    ode_lambda_functions += [sympy.utilities.lambdify(lambda_inputs, odes[i])]

input_vals = numpy.zeros(len(lambda_inputs))
input_vals[0] = decision_vector_values[0]
input_vals[1] = decision_vector_values[1]
input_vals[2] = decision_vector_values[2]
input_vals[3] = decision_vector_values[3]
input_vals[4] = decision_vector_values[4]
input_vals[5] = decision_vector_values[5]
input_vals[6] = decision_vector_values[6]
input_vals[7] = decision_vector_values[7]
input_vals[8] = decision_vector_values[8]

input_vals[13] = conservation_values[0]
input_vals[14] = conservation_values[1]
#input_vals[15] = 877876.0
#input_vals[15] = 878752.0
#input_vals[15] = conservation_values[2]

# solve the system dy/dt = f(y, t)
def f(y, t, inputs, ode_lambda_func):
    # s15 = y[0]
    # s3 = y[1]
    # s6 = y[2]
    # s16 = y[3]

    inputs[9] = y[0]#s15
    inputs[10] = y[1]#s3
    inputs[11] = y[2]#s6
    inputs[12] = y[3]#s16

    # the model equations
    f0 = ode_lambda_func[0](*tuple(inputs))
    f1 = ode_lambda_func[1](*tuple(inputs))
    f2 = ode_lambda_func[2](*tuple(inputs))
    f3 = ode_lambda_func[3](*tuple(inputs))

    return [f0, f1, f2, f3]

# time grid
t = numpy.linspace(1, 20., 10000)

#print(t)

# solve the ODEs
import scipy.integrate

# initial conditions
#y0 = [species_concentrations[6], species_concentrations[2], species_concentrations[3], species_concentrations[5]]
y0 = [34793.6, species_concentrations[2], species_concentrations[3], species_concentrations[5]]
input_vals[15] = 877876.0
#y0 = [60000.0, species_concentrations[2], species_concentrations[3], species_concentrations[5]]

soln = scipy.integrate.odeint(f, y0, t, args=(input_vals, ode_lambda_functions))

s15_sol1 = soln[:, 0]

y0 = [122611.0, species_concentrations[2], species_concentrations[3], species_concentrations[5]]
input_vals[15] = 878752.0
#y0 = [60000.0, species_concentrations[2], species_concentrations[3], species_concentrations[5]]

soln2 = scipy.integrate.odeint(f, y0, t, args=(input_vals, ode_lambda_functions))

s15_sol2 = soln2[:, 0]

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt


plt.plot(t, s15_sol1, label='Species s15 1')
plt.plot(t, s15_sol2, label='Species s15 2')
plt.ticklabel_format(axis='both', style='sci', scilimits=(-2, 2))
plt.show()
