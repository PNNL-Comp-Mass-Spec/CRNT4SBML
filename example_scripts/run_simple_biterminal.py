import sys
sys.path.insert(0, "..")
import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/simple_biterminal.xml")

network.print_c_graph()

approach = network.get_general_approach()

signal = "C2"
response = "s11"

approach.initialize_general_approach(signal=signal, response=response, fix_reactions=False)

bnds = [(2.4, 2.42), (27.5, 28.1), (2.0, 2.15), (48.25, 48.4), (0.5, 1.1), (1.8, 2.1), (17.0, 17.5), (92.4, 92.6),
        (0.01, 0.025), (0.2, 0.25), (0.78, 0.79), (3.6, 3.7), (0.15, 0.25), (0.06, 0.065)] + [(0.0, 100.0),
        (18.0, 18.5), (0.0, 100.0), (0.0, 100.0), (27.0, 27.1), (8.2, 8.3), (90.0, 90.1), (97.5, 97.9), (30.0, 30.1)]

# print(network.get_c_graph().get_species())
# print(approach.get_independent_species())
print(network.get_c_graph().get_ode_system())
print("")
print(approach.get_independent_odes())
# print(approach.get_independent_odes_subs())
# import sympy
# print("hi")
print(approach.get_conservation_laws())
#
#
# sys.exit()

iters = 15
d_iters = 1000
sd = 0
prnt_flg = False

# params_for_global_min, obj_fun_vals = approach.run_optimization(bounds=bnds, iterations=iters, seed=sd, print_flag=prnt_flg,
#                                                                 dual_annealing_iters=d_iters, confidence_level_flag=True, parallel_flag=True)
#
# if approach.get_my_rank() == 0:
#
#     import numpy
#     numpy.save('nuts_params.npy', params_for_global_min)

import numpy
params_for_global_min = numpy.load('nuts_params.npy')
print(len(approach.get_input_vector()))
print(approach.get_input_vector())
print(params_for_global_min)

# sys.exit()

list_of_ggplots = approach.run_direct_simulation(params_for_global_min, parallel_flag=True)

if approach.get_my_rank() == 0:

    print("list_of_ggplots ")
    g = list_of_ggplots[0]
    import plotnine as p9
    path = "./dir_sim_graphs"
    from matplotlib import rc
    rc('text', usetex=True)

    g = (g + p9.xlab("$E_{tot}$") + p9.ylab("$[S^{**}]$") + p9.scale_color_hue(labels=["High [$S^{**}$]", "Low [$S^{**}$]"]))
    g.save(filename=path + f"/sim_bif_diag_0_0.png", format="png", width=6, height=4, units='in', verbose=False)

    print("")

approach.generate_report()



