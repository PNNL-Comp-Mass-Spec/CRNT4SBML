########################################################################
########################################################################
# Please review the documentation provided at crnt4sbml.readthedocs.io #
# before running the code below.                                       #
########################################################################
########################################################################


import crnt4sbml
import numpy

network = crnt4sbml.CRNT("./sbml_files/Edelstein.xml")
signal = "C1"
response = "s5"
iters = 50

GA = network.get_general_approach()
GA.initialize_general_approach(signal=signal, response=response, fix_reactions=True)

bnds = [(1e-2, 1e2)]*len(GA.get_input_vector())

#######################################################################################################################
#######################################################################################################################
# The code below produces the values for bistability (in parallel): ####
########################################################################
# params_for_global_min, obj_fun_vals = GA.run_optimization(bounds=bnds, iterations=iters, seed=0, print_flag=False,
#                                                           dual_annealing_iters=1000, confidence_level_flag=True,
#                                                           parallel_flag=True)
# if GA.get_my_rank() == 0:
#     numpy.save('./optimization_parameters/edelstein_params.npy', params_for_global_min)
#######################################################################################################################
#######################################################################################################################


# loading in the parameters produced by optimization
params_for_global_min = numpy.load('./optimization_parameters/edelstein_params.npy')


#######################################################################################################################
#######################################################################################################################
# The code below produces Figure 10 of the paper (do not run in parallel!): ####
###############################################################################
multistable_param_ind, plot_specifications = GA.run_greedy_continuity_analysis(species=response,
                                                                               parameters=params_for_global_min,
                                                                               auto_parameters={'PrincipalContinuationParameter': signal},
                                                                               dir_path="./",
                                                                               plot_labels=["$B_{tot}$", "[A]", None])
#######################################################################################################################
#######################################################################################################################




#######################################################################################################################
#######################################################################################################################
# The code below produces Figures 11 and 12 of the paper Do not run in parallel!: ####
# Uncomment the below code and run it.                                            ####
######################################################################################

# def ff(t, cs, ks, ode_lambda_functions, jacobian):
#     return [i(*tuple(ks), *tuple(cs)) for i in ode_lambda_functions]
#
# def jac_f(t, cs, ks, ode_lambda_functions, jacobian):
#     return jacobian(*tuple(ks), *tuple(cs))
#
# param_index = 1
#
# C1_val = params_for_global_min[param_index][7] + params_for_global_min[param_index][8]
#
# y_fwd = [10.0, 0.0, C1_val]
# y_rev = [0.0, C1_val, 0.0]
#
# import scipy.integrate as itg
# import pandas
# from plotnine import ggplot, aes, geom_line, ylim, scale_color_distiller, facet_wrap, theme_bw, geom_path, \
#     geom_point, labs, annotate
# from plotnine import *
# from matplotlib import rc
# rc('text', usetex=True)
#
# kinetic_constants = params_for_global_min[param_index][0:6]
# end_t = 65000.0
# t = numpy.linspace(0.0, end_t, 32500)
# scan_vals = numpy.linspace(7.77, 7.805, 60)
#
# def sim_fun_fwd(x):
#     y_fwd[2] = x
#     return itg.solve_ivp(ff, (0.0, end_t), y_fwd, t_eval=t, method="BDF", args=(kinetic_constants, GA.get_ode_lambda_functions(),
#                                                        GA.get_jac_lambda_function()), jac=jac_f)
#
# def sim_fun_rev(x):
#     y_rev[1] = x
#     return itg.solve_ivp(ff, (0.0, end_t), y_rev, t_eval=t, method="BDF", args=(kinetic_constants, GA.get_ode_lambda_functions(),
#                                                        GA.get_jac_lambda_function()), jac=jac_f)
#
# sim_res_fwd = [numpy.transpose(sim_fun_fwd(i).y) for i in scan_vals]
# sim_res_rev = [numpy.transpose(sim_fun_rev(i).y) for i in numpy.flip(scan_vals)]
#
# out = pandas.DataFrame(columns=['dir', 'signal', 'time'] + network.get_c_graph().get_species())
#
# for i in range(len(sim_res_rev)):
#     out_i = pandas.DataFrame(sim_res_rev[i], columns=out.columns[3:])
#     out_i['time'] = t
#     out_i['signal'] = numpy.flip(scan_vals)[i]
#     out_i['dir'] = 'Low [A]'
#     out = pandas.concat([out, out_i[out.columns]])
#
# for i in range(len(sim_res_fwd)):
#     out_i = pandas.DataFrame(sim_res_fwd[i], columns=out.columns[3:])
#     out_i['time'] = t
#     out_i['signal'] = scan_vals[i]
#     out_i['dir'] = 'High [A]'
#     out = pandas.concat([out, out_i[out.columns]])
#
# out['dir'] = out['dir'].astype('category')
# out['dir'].cat.reorder_categories(['Low [A]', 'High [A]'], ordered=True, inplace=True)
#
# ###################### plotting ##################################
# g = (ggplot(out, aes('time', response, group='signal', color='signal'))
#      + geom_line(size=0.5)
#      + ylim(0, 5)
#      + labs(x="time", y="[A]")
#      + scale_color_distiller(palette='RdYlBu', type="diverging", name="$B_{tot}$")
#      + facet_wrap('~dir')
#      + theme_bw())
# g.save(filename="./Figure_11.png", format="png", width=8, height=4, units='in', verbose=False)
#
# out['dir'].cat.reorder_categories(['High [A]', 'Low [A]'], ordered=True, inplace=True)
#
# out = out[out.time == max(out.time)]
#
#
# def mask(df, key, value):
#     return df[df[key] == value]
#
#
# pandas.DataFrame.mask = mask
#
# # flip one of the directions
# xf = out.mask('dir', 'High [A]')
# xr = out.mask('dir', 'Low [A]')
#
# out = pandas.concat([xf, xr])
#
# g = (ggplot(out)
#          + aes(x='signal', y='s5', color='dir', alpha='dir', group='dir')
#          + geom_path(size=4, arrow=arrow(), show_legend=False)
#          + geom_point(size=2, shape='s', stroke=0.0)
#          + geom_point(color="black", size=2, alpha=1.0)
#          + scale_alpha_manual(values=[0.85, 0.35], guide=False)
#          + guides(color=guide_legend(override_aes={'size':6, 'alpha':[0.85, 0.35]}))
#          + theme(legend_title=element_blank(), text=element_text(size=14, weight='bold'),
#                  legend_key=element_rect(color='None', fill='None'),
#                  axis_text_x=element_line(color="black"),
#                  axis_text_y=element_line(color="black"),
#                  axis_title_x=element_text(size=19, weight='heavy'),
#                  axis_title_y=element_text(size=19, weight='heavy'),
#                  panel_background=element_rect(fill='None'), panel_border=element_rect(fill='None', color='#7f7f7f'),
#                  panel_grid_major=element_line(color='#E5E5E5', size=0.8),
#                  panel_grid_minor=element_line(color='#FAFAFA', size=1),
#                  strip_background=element_rect(fill='#CCCCCC', color='#7F7F7F', size=1)))
#
# g = (g + xlab("$B_{tot}$") + ylab("[A]")  + scale_color_manual(values=["red", "blue"], labels=["High [A]", "Low [A]"]))
#
# g.save(filename="./Figure_12.png", format="png", width=8, height=5, units='in', verbose=False)

#######################################################################################################################
#######################################################################################################################
