import sys
sys.path.insert(0, "..")
import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/Song.xml")

network.basic_report()

network.print_c_graph()

# network.get_network_graphml()

approach = network.get_general_approach()
approach.initialize_general_approach(signal="C1", response="s2", fix_reactions=True)

print(approach.get_input_vector())

bnds = [(1e-3, 6.0)]*len(network.get_c_graph().get_reactions()) + [(1e-3, 1000.0)]*len(network.get_c_graph().get_species())

# params, obj_fun_vals = approach.run_optimization(bounds=bnds, iterations=100, parallel_flag=True)
#
# if approach.get_my_rank() == 0:
#
#     import numpy
#     numpy.save('song_params.npy', params)

import numpy
params = numpy.load('song_params.npy')
params = [params[4]]

print(len(approach.get_input_vector()))
print(approach.get_input_vector())
print(params)

# sys.exit()

list_of_ggplots  = approach.run_direct_simulation(params, parallel_flag=True)

if approach.get_my_rank() == 0:

    print("list_of_ggplots ")
    g = list_of_ggplots[0]
    import plotnine as p9
    path = "./dir_sim_graphs"
    from matplotlib import rc
    rc('text', usetex=True)

    g = (g + p9.xlab("$E_{tot}$") + p9.ylab("$[S^*]$") + p9.scale_color_hue(labels=["High [$S^*$]", "Low [$S^*$]"]))
    g.save(filename=path + f"/sim_bif_diag_0_0.png", format="png", width=6, height=4, units='in', verbose=False)

    print("")

# multistable_param_ind, plot_specifications = approach.run_greedy_continuity_analysis(species="s2", parameters=params,
#                                                                                      auto_parameters={'PrincipalContinuationParameter': "C1"})

approach.generate_report()

