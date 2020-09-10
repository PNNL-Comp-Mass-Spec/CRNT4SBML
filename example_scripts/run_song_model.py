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

    # print("list_of_ggplots ")
    g = list_of_ggplots[0]
    import plotnine as p9
    path = "./dir_sim_graphs"
    from matplotlib import rc
    rc('text', usetex=True)

    # signal = '$E_{tot}$'
    # response = 's2'
    #
    # import pandas
    # from plotnine import *
    #
    # out = pandas.read_csv(path + '/song_dir_sim_data.txt', sep="\t", index_col=False)
    #
    # out = pandas.DataFrame(data=out, columns=['dir', 'signal', 's2'])
    #
    # print(f"out = {out}")
    #
    #
    # def mask(df, key, value):
    #     return df[df[key] == value]
    #
    #
    # pandas.DataFrame.mask = mask
    #
    # # flip one of the directions
    # xf = out.mask('dir', 'Forward scan') #out[out.dir].filter(like="Forward scan", axis=0)
    # xr = out.mask('dir', 'Reverse scan') #out.filter(like="Reverse scan", axis=0)
    #
    # print(f"xf = {xf}")
    #
    # print("")
    # print(f"xr = {xr}")
    #
    # xr = xr.iloc[::-1]
    #
    # print("")
    # print(f"xr = {xr}")
    #
    # # x2 = rbind(xf, xr[rev(seq_len(nrow(xr))),])
    #
    # out = pandas.concat([xf, xr])


    # sys.exit()

    # g = (ggplot(out)
    #      + aes(x='signal', y='s2', color='dir', alpha='dir', group='dir')
    #      + geom_path(size=4, arrow=arrow(), show_legend=False)
    #      + geom_point(size=2, shape='s', stroke=0.0)
    #      + geom_point(color="black", size=2, alpha=1.0)
    #      + scale_alpha_manual(values=[0.85, 0.35], guide=False)
    #      # + scale_color_manual(values=["red", "blue"], labels=["d", "k"])
    #      + guides(color=guide_legend(override_aes={'size':6, 'alpha':[0.85, 0.35]}))
    #      + theme(legend_title=element_blank(), text=element_text(size=14, weight='bold'),
    #              legend_key=element_rect(color='None', fill='None'),
    #              axis_text_x=element_line(color="black"),
    #              axis_text_y=element_line(color="black"),
    #              axis_title_x=element_text(size=19, weight='heavy'),
    #              axis_title_y=element_text(size=19, weight='heavy'),
    #              panel_background=element_rect(fill='None'), panel_border=element_rect(fill='None', color='#7f7f7f'),
    #              panel_grid_major=element_line(color='#E5E5E5', size=0.8),
    #              panel_grid_minor=element_line(color='#FAFAFA', size=1),
    #              strip_background=element_rect(fill='#CCCCCC', color='#7F7F7F', size=1)))

    g = (g + p9.xlab("$E_{tot}$") + p9.ylab("$[S^*]$")  + p9.scale_color_manual(values=["red", "blue"], labels=["High [$S^*$]", "Low [$S^*$]"]))


    g.save(filename=path + f"/song-bif-diag_2.png", format="png", width=8, height=5, units='in', verbose=False)


# multistable_param_ind, plot_specifications = approach.run_greedy_continuity_analysis(species="s2", parameters=params,
#                                                                                      auto_parameters={'PrincipalContinuationParameter': "C1"})

# approach.generate_report()

