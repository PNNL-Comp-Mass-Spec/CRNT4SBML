import sys
sys.path.insert(0, "..")
import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/Fig1Cii.xml")

approach = network.get_semi_diffusive_approach()

bounds = approach.get_optimization_bounds()

params_for_global_min, obj_fun_val_for_params = approach.run_optimization(bounds=bounds, parallel_flag=True)

approach.generate_report()
