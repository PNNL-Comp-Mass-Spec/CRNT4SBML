import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/Fig4B_open.xml")

approach = network.get_semi_diffusive_approach()

bounds = [(1e-2, 1e2)]*11

params_for_global_min, obj_fun_val_for_params = approach.run_optimization(bounds=bounds, iterations=10000,
                                                                          parallel_flag=True)

if approach.get_my_rank() == 0:

    network.basic_report()

    network.print_c_graph()

    ldt = network.get_low_deficiency_approach()
    ldt.report_deficiency_zero_theorem()
    ldt.report_deficiency_one_theorem()

    print("")

    approach.print_decision_vector()

    print("Key species:")
    print(approach.get_key_species())
    print("")

    print("Non key species:")
    print(approach.get_non_key_species())

    print("")
    print("Boundary species:")
    print(approach.get_boundary_species())

approach.generate_report()
