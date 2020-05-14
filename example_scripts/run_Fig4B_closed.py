import sys
sys.path.insert(0, "..")
import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/Fig4B_closed.xml")

approach = network.get_mass_conservation_approach()

bounds = [(1e-2, 1e2)]*11
concentration_bounds = [(1e-2, 1e2)]*3

params_for_global_min, obj_fun_val_for_params = approach.run_optimization(bounds=bounds,
                                                                          concentration_bounds=concentration_bounds,
                                                                          iterations=10000, parallel_flag=True)

if approach.get_my_rank() == 0:
    network.basic_report()

    network.print_c_graph()

    ldt = network.get_low_deficiency_approach()
    ldt.report_deficiency_zero_theorem()
    ldt.report_deficiency_one_theorem()

    print("")

    print("Decision Vector:")
    print(approach.get_decision_vector())
    print("")

    print("Species for concentration bounds:")
    print(approach.get_concentration_bounds_species())

    approach.generate_report()
