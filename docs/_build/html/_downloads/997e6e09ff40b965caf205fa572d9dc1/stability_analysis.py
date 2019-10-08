import dill
import sympy 
import numpy
import scipy.integrate
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
else:
    import matplotlib.pyplot as plt
import matplotlib.gridspec

#loading in values
with open("important_variables.dill", 'rb') as f:
    out = dill.load(f)

params_for_global_min = numpy.load('params.npy')

# Constructing necessary variables for ODE simulation
########################################################################################################################

# getting the ODE system
odes = out[0]

# construct sympy form of reactions and species
sympy_reactions = out[1]
sympy_species = out[2]

# creating the input for each ode lambda function
lambda_inputs = sympy_reactions + sympy_species

# creating a lambda function for each ODE to make evaluating them simpler
ode_lambda_functions = []
for i in range(len(odes)):
    ode_lambda_functions += [sympy.utilities.lambdify(lambda_inputs, odes[i])]

# first index of a species in input_vals
start_index = len(sympy_reactions)


# function to compute the derivative of y at t
def f(y, t, inputs, ode_lambda_func, start_ind):

    # setting species concentrations
    inputs[start_ind:] = y

    # the model equations
    ode_vals = []
    for i in ode_lambda_func:
        ode_vals.append(i(*tuple(inputs)))

    return ode_vals


########################################################################################################################


# assigning the values for ODE simulation
########################################################################################################################
########################################################################################################################

# setting the decision vector values that define the ODE simulation
decision_vector_values = params_for_global_min[2]

# choosing the plot specifications corresponding to the decision vector
plot_specs = out[6][0]

# setting time grid for ODE simulation
final_time = 150000  # index 2
# final_time = 2500 # index 3
t = numpy.linspace(0.0, final_time, 10000)

# setting initial species concentration
s15_init = 2.0e5 # index 2
# s15_init = 4.0e5 # index 3

# getting species concetrations for a specific parameter set
species_concentrations = []
for i in out[3]:
    species_concentrations.append(i(*tuple(decision_vector_values)))

# getting constant values corresponding to the conservation laws
conservation_values = out[4]*sympy.Matrix([species_concentrations]).T

# filling in the values for the reaction rates using optimization values
input_vals = numpy.zeros(len(lambda_inputs))
for i in range(len(sympy_reactions)):
    input_vals[i] = decision_vector_values[i]

########################################################################################################################
########################################################################################################################


# creating subplots to represent the ODE simulation graphically
########################################################################################################################
########################################################################################################################
########################################################################################################################


# Setting up the format of the subplot
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

# C3 values to vary in ode solution
C3_vec = numpy.linspace(plot_specs[0][0], plot_specs[0][1])

C3_min = min(C3_vec)
C3_max = max(C3_vec)

C1 = conservation_values[0]
C2 = conservation_values[1]

# constructing format for plots
fig = plt.figure(constrained_layout=True, figsize=(20, 6))
plt.rcParams.update({'font.size': 17})
spec = matplotlib.gridspec.GridSpec(ncols=3, nrows=1, figure=fig, width_ratios=[3.3, 3.3, 3.1], height_ratios=[1])

# setting up a colormap that's a simple transition
mymap = matplotlib.colors.LinearSegmentedColormap.from_list('mycolors', ['blue', 'red'])

# creating a colorbar that varies from C3_min to C3_max
CS3 = plt.cm.ScalarMappable(cmap=mymap, norm=plt.Normalize(vmin=C3_min, vmax=C3_max))

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

# Creating the middle subplot
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

fig.add_subplot(spec[0, 1])

steady_state_vals1 = numpy.zeros(len(C3_vec))
for i in range(len(C3_vec)):
    y0 = [C3_vec[i], C2, 0.0, 0.0, C1, 0.0, 0.0]
    soln1 = scipy.integrate.odeint(f, y0, t, args=(input_vals, ode_lambda_functions, start_index))
    s15_sol = soln1[:, 6]
    steady_state_vals1[i] = soln1[-1, 6]
    r = (C3_vec[i] - C3_min)/(C3_max - C3_min)
    plt.plot(t, s15_sol, color=(r, 0, 1-r))

plt.xlabel("time (seconds)")
plt.ylim(plot_specs[1])
plt.ticklabel_format(axis='both', style='sci', scilimits=(-2, 2))
plt.title("Initial [s15] = 0 pM")

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

# Creating the rightmost subplot
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

fig.add_subplot(spec[0, 2])

steady_state_vals2 = numpy.zeros(len(C3_vec))
for i in range(len(C3_vec)):
    y0 = [C3_vec[i]-2*s15_init, C2, 0.0, 0.0, C1, 0.0, s15_init]
    soln1 = scipy.integrate.odeint(f, y0, t, args=(input_vals, ode_lambda_functions, start_index))
    s15_sol = soln1[:, 6]
    steady_state_vals2[i] = soln1[-1, 6]
    r = (C3_vec[i] - C3_min)/(C3_max - C3_min)
    plt.plot(t, s15_sol, color=(r, 0, 1-r))

plt.xlabel("time (seconds)")
plt.ylim(plot_specs[1])
plt.ticklabel_format(axis='both', style='sci', scilimits=(-2, 2))
plt.title("Initial [s15] = " + numpy.format_float_scientific(s15_init, exp_digits=0, trim='-').replace('+', '') + " pM")

clb = plt.colorbar(CS3)
clb.set_label("C3")
clb.formatter.set_powerlimits((2, 2))
clb.ax.ticklabel_format(axis='both', style='sci', scilimits=(-2, 2))
clb.update_ticks()

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

# Creating the leftmost subplot
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

fig.add_subplot(spec[0, 0])

plt.plot(C3_vec, steady_state_vals1, 'bo')
plt.plot(C3_vec, steady_state_vals2, 'bo')
plt.xlabel("C3")
plt.ylabel("[s15] pM")
plt.xlim(plot_specs[0])
plt.ylim(plot_specs[1])
plt.ticklabel_format(axis='both', style='sci', scilimits=(-2, 2))
plt.title("[s15] at " + numpy.format_float_scientific(final_time, exp_digits=0, trim='-').replace('+', '') + " seconds")


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

plt.savefig('bistability_of_ODEs.png')
plt.close(fig)

########################################################################################################################
########################################################################################################################
########################################################################################################################
