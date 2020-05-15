import sys
sys.path.insert(0, "..")
import crnt4sbml

network = crnt4sbml.CRNT("../sbml_files/Fig1Ci.xml")

approach = network.get_mass_conservation_approach()

bounds, concentration_bounds = approach.get_optimization_bounds()

params_for_global_min, obj_fun_val_for_params = approach.run_optimization(bounds=bounds, concentration_bounds=concentration_bounds)

approach.run_direct_simulation(response="s15", signal="C3", params_for_global_min=params_for_global_min)

approach.generate_report()
