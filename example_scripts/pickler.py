import crnt4sbml_test
import pickle
import contextlib
import io

def generate_c_graph_list(crnt_obj, flag):
    a = crnt_obj.get_c_graph().get_a()
    complexes = crnt_obj.get_c_graph().get_complexes()
    deficiency = crnt_obj.get_c_graph().get_deficiency()
    dim_eq_man = crnt_obj.get_c_graph().get_dim_equilibrium_manifold()
    g_edges = crnt_obj.get_c_graph().get_g_edges()
    g_nodes = crnt_obj.get_c_graph().get_g_nodes()
    reversible = crnt_obj.get_c_graph().get_if_cgraph_weakly_reversible()
    link_def = crnt_obj.get_c_graph().get_linkage_classes_deficiencies()
    classification = crnt_obj.get_c_graph().get_network_dimensionality_classification()
    num_term = crnt_obj.get_c_graph().get_number_of_terminal_strong_lc_per_lc()
    ode_system = crnt_obj.get_c_graph().get_ode_system()
    psi = crnt_obj.get_c_graph().get_psi()
    reactions = crnt_obj.get_c_graph().get_reactions()
    s = crnt_obj.get_c_graph().get_s()
    species = crnt_obj.get_c_graph().get_species()
    revers_link = crnt_obj.get_c_graph().get_weak_reversibility_of_linkage_classes()
    y = crnt_obj.get_c_graph().get_y()

    if flag:
        b = crnt_obj.get_c_graph().get_b()
        lambda_val = crnt_obj.get_c_graph().get_lambda()

        return [a, b, complexes, deficiency, dim_eq_man, g_edges, g_nodes, reversible, link_def, classification,
                num_term, ode_system, psi, reactions, s, species, revers_link, y, lambda_val]
    else:
        return [a, complexes, deficiency, dim_eq_man, g_edges, g_nodes, reversible, link_def, classification,
                num_term, ode_system, psi, reactions, s, species, revers_link, y]


def generate_low_def_list(crnt_obj):

    d_s_any = crnt_obj.get_low_deficiency_approach().does_satisfy_any_low_deficiency_theorem()
    d_s_one = crnt_obj.get_low_deficiency_approach().does_satisfy_deficiency_one_theorem()
    d_s_zero = crnt_obj.get_low_deficiency_approach().does_satisfy_deficiency_zero_theorem()

    f = io.StringIO()
    with contextlib.redirect_stdout(f):
        crnt_obj.get_low_deficiency_approach().report_deficiency_one_theorem()
    report_one = f.getvalue()
    f.close()

    f = io.StringIO()
    with contextlib.redirect_stdout(f):
        crnt_obj.get_low_deficiency_approach().report_deficiency_zero_theorem()
    report_zero = f.getvalue()
    f.close()

    return [d_s_any, d_s_one, d_s_zero, report_one, report_zero]

def generate_mass_con_list(crnt_obj):

    approach = crnt_obj.get_mass_conservation_approach()
    concentration_solutions = approach.get_concentration_solutions()
    concentration_vals = approach.get_concentration_vals()
    conservation_laws = approach.get_conservation_laws()
    dch_matrix = approach.get_dch_matrix()
    decision_vector = approach.get_decision_vector_x()
    g_matrix = approach.get_g_matrix()
    h_vector = approach.get_h_vector()
    objective_fun_params = approach.get_objective_fun_params()
    symbolic_objective_fun = approach.get_symbolic_objective_fun()
    w_matrix = approach.get_w_matrix()
    w_nullspace = approach.get_w_nullspace()

    return [concentration_solutions, concentration_vals, conservation_laws, dch_matrix, decision_vector, g_matrix,
            h_vector, objective_fun_params, symbolic_objective_fun, w_matrix, w_nullspace]

def generate_semi_diff_list(crnt_obj):

    approach = crnt_obj.get_semi_diffusive_approach()
    boundary_species = approach.get_boundary_species()
    decision_vector = approach.get_decision_vector()
    key_species = approach.get_key_species()
    mu_vector = approach.get_mu_vector()
    non_key_species = approach.get_non_key_species()
    s_to_matrix = approach.get_s_to_matrix()
    symbolic_objective_fun = approach.get_symbolic_objective_fun()
    symbolic_polynomial_fun = approach.get_symbolic_polynomial_fun()
    y_r_matrix = approach.get_y_r_matrix()

    f = io.StringIO()
    with contextlib.redirect_stdout(f):
        approach.print_decision_vector()
    print_decision_vector = f.getvalue()
    f.close()

    return [boundary_species, decision_vector, key_species, non_key_species, mu_vector, s_to_matrix,
            symbolic_objective_fun, symbolic_polynomial_fun, y_r_matrix, print_decision_vector]


#mass_con_crnt_obj = crnt4sbml_test.CRNT("../sbml_files/Fig1Ci.xml")
#mass_con_c_graph_list = generate_c_graph_list(mass_con_crnt_obj, True)
#mass_con_low_def_list = generate_low_def_list(mass_con_crnt_obj)
#mass_con_list = generate_mass_con_list(mass_con_crnt_obj)

semi_diff_crnt_obj = crnt4sbml_test.CRNT("../sbml_files/Fig1Cii.xml")
#semi_diff_c_graph_list = generate_c_graph_list(semi_diff_crnt_obj, False)
#semi_diff_low_def_list = generate_low_def_list(semi_diff_crnt_obj)
semi_diff_list = generate_semi_diff_list(semi_diff_crnt_obj)
#print(semi_diff_list)


#low_def_crnt_obj = crnt4sbml_test.CRNT("../sbml_files/feinberg_ex3_13.xml")
#low_def_c_graph_list = generate_c_graph_list(low_def_crnt_obj, False)
#low_def_low_def_list = generate_low_def_list(low_def_crnt_obj)


#with open('mass_con_list.pickle', 'wb') as outf:
#    outf.write(pickle.dumps(mass_con_list))

with open('semi_diff_list.pickle', 'wb') as outf:
    outf.write(pickle.dumps(semi_diff_list))

#with open('low_def_low_def_list.pickle', 'wb') as outf:
#    outf.write(pickle.dumps(low_def_low_def_list))

