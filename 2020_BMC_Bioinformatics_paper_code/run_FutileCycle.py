########################################################################
########################################################################
# Please review the documentation provided at crnt4sbml.readthedocs.io #
# before running the code below.                                       #
########################################################################
########################################################################


import crnt4sbml
import numpy

network = crnt4sbml.CRNT("./sbml_files/FutileCycle.xml")

approach = network.get_general_approach()
approach.initialize_general_approach(signal="C1", response="s2", fix_reactions=True)

bnds = [(1e-3, 6.0)]*len(network.get_c_graph().get_reactions()) + \
       [(1e-3, 1000.0)]*len(network.get_c_graph().get_species())

#######################################################################################################################
#######################################################################################################################
# The code below produces the values for bistability (in parallel): ####
########################################################################
# params, obj_fun_vals = approach.run_optimization(bounds=bnds, iterations=100, parallel_flag=True)
#
# if approach.get_my_rank() == 0:
#
#     import numpy
#     numpy.save('./optimization_parameters/FutileCycle.npy', params)
#######################################################################################################################
#######################################################################################################################


# loading in the parameters produced by optimization
params = numpy.load('./optimization_parameters/FutileCycle.npy')
params = [params[4]] # specific parameter values used in the paper


#######################################################################################################################
#######################################################################################################################
# The code below produces Figure 5 (in parallel): ####
######################################################
list_of_ggplots  = approach.run_direct_simulation(params, parallel_flag=True)

if approach.get_my_rank() == 0:

    g = list_of_ggplots[0]
    import plotnine as p9
    from matplotlib import rc
    rc('text', usetex=True)

    g = (g + p9.xlab("$E_{tot}$") + p9.ylab("$[S^*]$")  + p9.scale_color_manual(values=["red", "blue"],
                                                                                labels=["High [$S^*$]",
                                                                                        "Low [$S^*$]"]))

    g.save(filename=f"./Figure_5.png", format="png", width=8, height=5, units='in', verbose=False)

#######################################################################################################################
#######################################################################################################################
