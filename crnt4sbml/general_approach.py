import sympy
import sys
import itertools
import numpy
import scipy.optimize
import scipy.linalg
import warnings
import time
import os
from .bistability_finder import BistabilityFinder


class GeneralApproach:
    """
    Class for constructing a more general approach to bistability detection for systems with mass action kinetics.
    """
    def __init__(self, cgraph, signal, response):
        """
        Initialization of GeneralApproach class using a network contructed by crnt4sbml.CRNT().

        See also
        ---------
        fill in this
        """

        self.__cgraph = cgraph

        self.__signal = signal
        self.__response = response

        if not self.__cgraph.get_dim_equilibrium_manifold() > 0:
            print("# of species - rank(S) is not greater than zero!")
            print("The general approach cannot be ran!")
            sys.exit()

        # todo: check that there is more than one independent ODE

        self.__full_system = self.__cgraph.get_ode_system()

        self.__sympy_species = [sympy.Symbol(i, positive=True) for i in self.__cgraph.get_species()]
        self.__sympy_reactions = [sympy.Symbol(i, positive=True) for i in self.__cgraph.get_reactions()]

        self.__N = len(self.__sympy_species)
        self.__R = len(self.__sympy_reactions)
        self.__M = len(self.__cgraph.get_complexes())
        self.__BT_mat = None

        self.__A_mat = self.__cgraph.get_a()
        self.__Y_mat = self.__cgraph.get_y()
        self.__psi_vec = self.__cgraph.get_psi()
        self.__S_mat = self.__cgraph.get_s()

        try:
            self.__BT_mat = self.__cgraph.get_b()
        except:
            self.__create_BT_matrix()
            sympy.pprint(self.__BT_mat)

        self.__cons_laws_sympy = None
        self.__indp_species = None
        self.__indp_system = None
        self.__indp_system_subs = None
        self.__replacements = None
        self.__det_jac_indp_subs = None
        self.__D2_g_u_w = None
        self.__decision_vector = None
        self.__lambda_indp_odes = None
        self.__lambda_det_jac = None
        self.__lambda_jec = None
        self.__lambda_D2_g_u_w = None
        self.__lambda_jac_eps = None
        self.__important_info = ""
        self.__lambda_variables = None
        self.__indices_in_dec_vec = None

        self.__create_conservation_laws()
        self.__construct_independent_system()
        self.__construct_independent_system_with_substitution()
        self.__create_variables_for_optimization()

    def __create_BT_matrix(self):
        # (bt = NullSpace[Transpose[Y.A]]) // MatrixForm
        the_null_space = (self.__Y_mat * self.__A_mat).T.nullspace()

        # getting the number of columns given in nullspace computation
        sizes = len(the_null_space)

        # simplifying the entries given in the columns
        # produced by the nullspace
        for i in range(0, sizes):
            the_null_space[i] = sympy.simplify(the_null_space[i])

        bt_temp = self.__create_non_negative_b_matrix(the_null_space, sizes)

        # taking out zero rows if they are in unique vectors
        bt = sympy.Matrix(bt_temp[~numpy.all(bt_temp == 0, axis=1)]).evalf(10)

        # making the dimension of bt be lambda by N, just incase
        # we have more vectors than we should
        bt = bt[0:sizes, :]

        # putting in a check to make sure that the rows
        # of bt are linearly independent
        _, vals = bt.T.rref()
        if len(vals) != sizes:
            print(" ")
            print("Error: Not all rows of B^T are linearly independent!")
            print(" ")
            sys.exit()

        self.__BT_mat = bt

    def __create_non_negative_b_matrix(self, the_null_space, sizes):

        a_null = numpy.zeros((self.__N, sizes))
        for i in range(self.__N):
            temp_vec = []
            for j in range(sizes):
                temp_vec.append(the_null_space[j][i])
            a_null[i, :] = temp_vec

        a_eq = numpy.array([numpy.sum(a_null, axis=0)])
        # must multiply by negative one because in optimization we have the inequality <= 0.0
        a_ub = -1.0 * a_null
        b_ub = numpy.zeros(self.__N)
        b_eq = numpy.array([1.0])

        # defining the number of solutions to simulate
        num_sols = ((a_ub.shape[0] + 1) * (a_ub.shape[1] + 1)) * 10

        # a matrix to hold the different solutions of the method
        sol = numpy.zeros((num_sols, self.__N))

        # searching the nullspace for nonnegative vectors
        for i in range(num_sols):
            # generating a multivariate normal random variates
            crit = numpy.random.normal(0, 1.0, a_ub.shape[1])

            # Solving the linear programming problem
            out = scipy.optimize.linprog(crit, A_eq=a_eq, b_eq=b_eq, A_ub=a_ub, b_ub=b_ub, bounds=(None, None),
                                         method='simplex')
            # multiplying our nullspace vectors with the minimized value
            # to create a nonnegative vector
            sol[i, :] = numpy.dot(a_null, out.x)
            sol[i, [numpy.abs(sol[i, :]) < numpy.finfo(float).eps][0]] = numpy.float64(0.0)

        # getting the unique vectors
        unique_vecs = numpy.unique(sol, axis=0)
        num_rows = numpy.size(unique_vecs, 0)

        # taking the smallest nonzero entry of the unique vectors
        # and dividing by it this will hopfully create nice
        # looking vectors that are whole numbers
        for i in range(num_rows):
            minval = numpy.min(unique_vecs[i, numpy.nonzero(unique_vecs[i, :])])
            unique_vecs[i, :] = unique_vecs[i, :] / numpy.float64(minval)

        unique_vecs = numpy.unique(unique_vecs.round(decimals=10), axis=0)

        return unique_vecs

    def __create_conservation_laws(self):

        self.__cons_laws_sympy = self.__BT_mat*sympy.Matrix(self.__sympy_species)
        # print("conservation laws")
        # sympy.pprint(self.__cons_laws_sympy)

    def __construct_independent_system(self):

        # forming ya matrix
        ya = self.__Y_mat * self.__A_mat

        # finding how many rows are indep in ya
        _, vals = ya.T.rref()

        num_indp_eqns = len(vals)
        num_dep_eqns = self.__N - num_indp_eqns

        # getting dimensions of bt
        bt_rows = self.__BT_mat.shape[0]
        bt_cols = self.__BT_mat.shape[1]

        bt_nonzero_ind = []
        species_num = self.__sympy_species.index(sympy.Symbol(self.__response, positive=True))
        for i in range(bt_rows):
            bt_nonzero_ind.append([j for j in range(bt_cols) if self.__BT_mat[i, j] != 0 and j != species_num])

        chosen_indp_indices, chosen_dep_indices = self.__get_indp_dep_species_indices(bt_nonzero_ind, num_dep_eqns,
                                                                                      num_indp_eqns, ya)

        self.__replacements, self.__indp_species, self.__indp_system = self.__construct_important_variables(chosen_indp_indices,
                                                                                                            chosen_dep_indices, ya)

    def __construct_independent_system_with_substitution(self):

        self.__indp_system_subs = sympy.zeros(len(self.__indp_system), 1)
        self.__indp_system_subs = self.__indp_system[:, :]

        # making the replacements in the indep. ODEs
        for i in range(len(self.__indp_system_subs)):
            for j in range(len(self.__replacements)):
                self.__indp_system_subs[i] = self.__indp_system_subs[i].subs(self.__replacements[j][0], self.__replacements[j][1])

    def __get_indp_dep_species_indices(self, bt_nonzero_ind, num_dep_eqns, num_indp_eqns, ya):
        # getting all combinations of the list indices
        possible_dep_species = list(itertools.product(*bt_nonzero_ind))

        removed_entries = []
        # remove tuples that have duplicate entries
        for i in range(len(possible_dep_species)):
            if len(set(possible_dep_species[i])) != num_dep_eqns:
                removed_entries.append(i)
        for index in sorted(removed_entries, reverse=True):
            del possible_dep_species[index]

        # get corresponding possible dependent species
        possible_indp_species = []
        species_ind = [i for i in range(self.__N)]

        for i in possible_dep_species:
            possible_indp_species.append([j for j in species_ind if j not in i])

        # using YA to pick one of the possible indices
        chosen_indp_indices = []
        chosen_dep_indices = []
        for i in range(len(possible_indp_species)):
            _, vals = ya[possible_indp_species[i], :].T.rref()

            if len(vals) == num_indp_eqns:
                chosen_indp_indices = possible_indp_species[i]
                chosen_dep_indices = possible_dep_species[i]
                break

        return chosen_indp_indices, chosen_dep_indices

    def __construct_important_variables(self, chosen_indp_indices, chosen_dep_indices, ya):
        # getting independent concentrations
        ind_spec_conc_temp = [self.__sympy_species[i] for i in chosen_indp_indices]

        # getting dependent concentrations
        dep_spec_conc = [self.__sympy_species[i] for i in chosen_dep_indices]

        # constructing the independent ODEs
        indp_odes_temp = ya[chosen_indp_indices, :] * self.__psi_vec

        # Lambda function of conservation laws
        self.__cons_laws_lamb = [sympy.utilities.lambdify(self.__sympy_species, self.__cons_laws_sympy[i])
                                 for i in range(len(self.__cons_laws_sympy))]

        cons_laws_sympy_eq = [sympy.Eq(sympy.Symbol('C' + str(i + 1), real=True), self.__cons_laws_sympy[i])
                              for i in range(len(self.__cons_laws_sympy))]

        dep_conc_in_laws = self.__dependent_species_concentrations(self.__cons_laws_sympy, dep_spec_conc)

        replacements = self.__find_dep_concentration_replacements(dep_conc_in_laws, self.__cons_laws_sympy,
                                                                  dep_spec_conc, cons_laws_sympy_eq)

        return replacements, ind_spec_conc_temp, indp_odes_temp

    def __dependent_species_concentrations(self, cons_laws_sympy, dep_spec_conc):
        # finding those dep_spec_conc that occur in each conservation law
        dep_conc_in_laws = []
        for i in range(len(cons_laws_sympy)):
            temp = []
            for j in range(len(dep_spec_conc)):
                if cons_laws_sympy[i].count(dep_spec_conc[j]) > 0:
                    temp.append(dep_spec_conc[j])
            dep_conc_in_laws.append(temp)
        return dep_conc_in_laws

    def __is_list_empty(self, inlist):
        if isinstance(inlist, list):  # Is a list
            return all(map(self.__is_list_empty, inlist))
        return False  # Not a list

    def __find_dep_concentration_replacements(self, dep_conc_in_laws, cons_laws_sympy, dep_spec_conc,
                                              cons_laws_sympy_eq):
        replacements = []
        flag = True
        while flag:
            for i in range(len(cons_laws_sympy_eq)):
                if len(dep_conc_in_laws[i]) == 1:
                    temp = sympy.solve(cons_laws_sympy_eq[i], dep_conc_in_laws[i])

                    cons_laws_sympy = [cons_laws_sympy[j].subs(dep_conc_in_laws[i][0], temp[0])
                                       for j in range(len(cons_laws_sympy))]

                    cons_laws_sympy_eq = [sympy.Eq(sympy.Symbol('C' + str(i+1), real=True), cons_laws_sympy[i])
                                          for i in range(len(cons_laws_sympy))]

                    replacements.append([dep_conc_in_laws[i][0], temp[0]])

            dep_conc_in_laws = self.__dependent_species_concentrations(cons_laws_sympy, dep_spec_conc)

            if self.__is_list_empty(dep_conc_in_laws):
                flag = False

        return replacements

    def __create_variables_for_optimization(self):

        self.__jac_indp_subs = self.__indp_system_subs.jacobian(sympy.Matrix(self.__indp_species))

        u_variables = sympy.Matrix(self.__indp_species)
        w_vector = sympy.Matrix([sympy.Symbol('w' + str(i)) for i in range(len(u_variables))])

        self.__D2_g_u_w = self.second_derivative_operator(self.__indp_system_subs, u_variables, w_vector)

        self.__det_jac_indp_subs = self.__jac_indp_subs.det(method='lu')

        self.__lambda_variables = self.__sympy_reactions + self.__indp_species + \
                                 [sympy.Symbol("C" + str(i+1), real=True) for i in range(len(self.__cons_laws_sympy))]

        jac_eps = self.__indp_system_subs.jacobian([sympy.Symbol(self.__signal, real=True)])

        self.__lambda_jac_eps = sympy.utilities.lambdify(self.__lambda_variables, jac_eps)

        self.__lambda_indp_odes = [sympy.utilities.lambdify(self.__lambda_variables, i) for i in self.__indp_system_subs]

        self.__lambda_det_jac = sympy.utilities.lambdify(self.__lambda_variables, self.__det_jac_indp_subs)

        self.__lambda_jec = sympy.utilities.lambdify(self.__lambda_variables, self.__jac_indp_subs)

        self.__lambda_D2_g_u_w = sympy.utilities.lambdify(self.__lambda_variables + [sympy.Symbol('w' + str(i))
                                                                                    for i in range(len(u_variables))],
                                                          self.__D2_g_u_w)

        self.__enforce_steady_state()

        self.__decision_vector = [i for i in self.__sympy_reactions if i not in self.__fixed_reactions] + self.__sympy_species

    @staticmethod
    def second_derivative_operator(g_functions, u_variables, w_vector):

        Dg_u = g_functions.jacobian(u_variables)

        # operating on w_vector
        Dg_u_w = Dg_u*w_vector

        D2_g_u_w = Dg_u_w.jacobian(u_variables)

        return D2_g_u_w

    @staticmethod
    def unique_everseen(iterable, key=None):
        "List unique elements, preserving order. Remember all elements ever seen."
        # unique_everseen('AAAABBBCCDAABBB') --> A B C D
        # unique_everseen('ABBCcAD', str.lower) --> A B C D
        seen = set()
        seen_add = seen.add
        if key is None:
            for element in itertools.filterfalse(seen.__contains__, iterable):
                seen_add(element)
                yield list(element)
        else:
            for element in iterable:
                k = key(element)
                if k not in seen:
                    seen_add(k)
                    yield list(element)

    def __construct_possible_fixed_reactions(self):

        # obtaining the reactions in each ODE
        reactions_lists = []
        for i in self.__indp_system_subs:
            atoms = list(i.atoms(sympy.Symbol))
            reactions = [j for j in self.__sympy_reactions if j in atoms]
            reactions_lists.append(reactions)

        # finding all combinations of the reactions for each ODE minus those with repeating elements
        temp = [p for p in itertools.product(*reactions_lists) if len(set(p)) == len(p)]

        # finds unique elements of the list where the elements are seen as unordered sets
        possible_fixed_reactions = list(self.unique_everseen(temp, key=frozenset))

        return possible_fixed_reactions

    def __enforce_steady_state(self):

        self.__fixed_reactions = []
        a, vals = self.__S_mat.rref()

        self.__fixed_reactions = [self.__sympy_reactions[i] for i in vals]

        if self.__fixed_reactions:
            temp_solution_tuple = sympy.solve(self.__indp_system_subs, self.__fixed_reactions, dict=True)
            if len(temp_solution_tuple) == 1:
                self.__soln_to_fixed_reactions = [sympy.factor(temp_solution_tuple[0].get(i)) for i in self.__fixed_reactions]
            else:
                # todo: check for positive solutions here
                print("multiple solutions found in reaction solve.")
                self.__soln_to_fixed_reactions = []

        self.__vars_for_lam_fixed_reactions = [i for i in self.__lambda_variables if i not in self.__fixed_reactions]
        self.__lambda_fixed_reactions = [sympy.utilities.lambdify(self.__vars_for_lam_fixed_reactions, i) for i in self.__soln_to_fixed_reactions]

    def __feasible_point_obj_func(self, x, boundz, full_set_of_values, reaction_ind, cons_laws_lamb, indp_spec_ind, inputs):

        # filling in values for non-fixed reactions
        inputs[0:self.__R - len(reaction_ind)] = x[0:self.__R - len(reaction_ind)]

        # filling in values for independent species
        count = self.__R - len(reaction_ind)
        for i in indp_spec_ind:
            inputs[count] = x[i]
            count += 1

        # filling in values for conservation laws
        for i in cons_laws_lamb:
            inputs[count] = i(*tuple(x[self.__R - len(reaction_ind):]))
            count += 1

        # filling in all species values
        full_set_of_values[self.__R:] = x[self.__R - len(reaction_ind):]

        # filling in non-fixed reaction values
        count = 0
        for i in range(self.__R):
            if i not in reaction_ind:
                full_set_of_values[i] = x[count]
                count += 1

        # filling in fixed reaction values
        full_set_of_values[reaction_ind] = numpy.array([i(*tuple(inputs)) for i in self.__lambda_fixed_reactions])

        sumval = numpy.float64(0.0)
        for i in range(len(full_set_of_values)):
            sumval += numpy.maximum(numpy.float64(0.0), numpy.float64(boundz[i][0]) - full_set_of_values[i]) ** 2
            sumval += numpy.maximum(numpy.float64(0.0), full_set_of_values[i] - numpy.float64(boundz[i][1])) ** 2

        return sumval

    def __obj_func(self, x, jac_det, boundz, full_set_of_values, reaction_ind, cons_laws_lamb, indp_spec_ind, inputs, input_values):

        # filling in values for non-fixed reactions
        inputs[0:self.__R - len(reaction_ind)] = x[0:self.__R - len(reaction_ind)]

        # filling in values for independent species
        count = self.__R - len(reaction_ind)
        for i in indp_spec_ind:
            inputs[count] = x[i]
            count += 1

        # filling in values for conservation laws
        for i in cons_laws_lamb:
            inputs[count] = i(*tuple(x[self.__R - len(reaction_ind):]))
            count += 1

        # filling in all species values
        full_set_of_values[self.__R:] = x[self.__R - len(reaction_ind):]

        # filling in non-fixed reaction values
        count = 0
        for i in range(self.__R):
            if i not in reaction_ind:
                full_set_of_values[i] = x[count]
                count += 1

        # filling in fixed reaction values
        full_set_of_values[reaction_ind] = numpy.array([i(*tuple(inputs)) for i in self.__lambda_fixed_reactions])

        # filing in reactions
        input_values[0: self.__R] = full_set_of_values[0:self.__R]

        # filling in independent species
        input_values[self.__R: self.__R + len(indp_spec_ind)] = x[indp_spec_ind]

        # filling in conservation laws
        input_values[self.__R + len(indp_spec_ind):] = inputs[self.__R - len(reaction_ind)+ len(indp_spec_ind):]

        fun = numpy.float64(0.0)

        determinant = (jac_det(*tuple(input_values)))**2

        if numpy.isfinite(determinant):
            fun += determinant

        else:
            fun += numpy.PINF

        sumval = numpy.float64(0.0)
        for i in range(len(full_set_of_values)):
            sumval += numpy.maximum(numpy.float64(0.0), numpy.float64(boundz[i][0]) - full_set_of_values[i]) #** 2
            sumval += numpy.maximum(numpy.float64(0.0), full_set_of_values[i] - numpy.float64(boundz[i][1])) #** 2

        fun += sumval

        return fun

    def __feasible_point_method(self, bounds, num_constraint_method_iters, seed, print_flag):

        all_vars_with_bounds = self.__sympy_reactions + self.__sympy_species
        decision_vector_bounds_ind = [all_vars_with_bounds.index(i) for i in self.__decision_vector]

        decision_vector_bounds = [bounds[i] for i in decision_vector_bounds_ind]

        samples = numpy.random.rand(num_constraint_method_iters, len(self.__decision_vector))

        ranges = numpy.asarray(decision_vector_bounds, dtype=numpy.float64)
        samples = samples * (ranges[:, 1] - ranges[:, 0]) + ranges[:, 0]

        full_set_of_values = numpy.zeros(len(bounds), dtype=numpy.float64)

        fixed_reaction_ind_all = [all_vars_with_bounds.index(i) for i in self.__fixed_reactions]
        indp_spec_ind_dec = [self.__decision_vector.index(i) for i in self.__indp_species]

        inputs = numpy.zeros(len(self.__vars_for_lam_fixed_reactions), dtype=numpy.float64)

        feasible_point_sets = []
        for i in range(num_constraint_method_iters):

            with numpy.errstate(divide='ignore', invalid='ignore'):

                # result = scipy.optimize.basinhopping(self.__feasible_point_obj_func, samples[i],
                #                                      minimizer_kwargs={'method': 'Nelder-Mead',
                #                                                        'args': (bounds, full_set_of_values,
                #                                                                 fixed_reaction_ind_all,
                #                                                                 self.__cons_laws_lamb,
                #                                                                 indp_spec_ind_dec, inputs),
                #                                                        'tol': 1e-16},
                #                                      niter=2, seed=seed)

                result = scipy.optimize.dual_annealing(self.__feasible_point_obj_func, bounds=decision_vector_bounds,
                                                       args=(bounds, full_set_of_values, fixed_reaction_ind_all,
                                                             self.__cons_laws_lamb, indp_spec_ind_dec, inputs),
                                                       x0=samples[i],
                                                       seed=seed,
                                                       local_search_options={'method': "Nelder-Mead"}, maxiter=100)

                if print_flag:
                    print("Objective function value: " + str(result.fun))
                    print("Decision vector used: ")
                    print(result.x)

                if abs(result.fun) > numpy.float64(1e-100):
                    result1 = scipy.optimize.minimize(self.__feasible_point_obj_func, result.x,
                                                      args=(bounds, full_set_of_values, fixed_reaction_ind_all,
                                                            self.__cons_laws_lamb, indp_spec_ind_dec, inputs),
                                                      method='Nelder-Mead', tol=1e-16)

                    if print_flag:
                        print("Objective function value: " + str(result1.fun))
                        print("Decision vector used: ")
                        print(result1.x)

                    if abs(result1.fun) <= numpy.finfo(float).eps:
                        feasible_point_sets.append(result1.x)

                else:
                    feasible_point_sets.append(result.x)

                if print_flag:
                    print("")

        return feasible_point_sets, fixed_reaction_ind_all, indp_spec_ind_dec, decision_vector_bounds

    def run_optimization(self, bounds=None, iterations=10, seed=0, print_flag=False, dual_annealing_iters=1000):

        """
                Function for running the optimization problem for the mass conservation approach.

                Parameters
                -----------
                    bounds: list of tuples
                        A list defining the lower and upper bounds for each variable in the decision vector.
                    iterations: int
                        The number of iterations to run the feasible point method.
                    seed: int
                        Seed for the random number generator. None should be used if a random generation is desired.
                    print_flag: bool
                        Should be set to True if the user wants the objective function values found in the optimization problem
                        and False otherwise.
                Returns
                --------
                params_for_global_min: list of numpy arrays
                    A list of numpy arrays that correspond to the decision vectors of the problem.
                obj_fun_val_for_params: list of floats
                    A list of objective function values produced by the corresponding decision vectors in params_for_global_min.

                Examples
                ---------
                Fill in
                """

        print("Starting optimization ...")
        start_t = time.time()
        numpy.random.seed(seed)

        print("Starting feasible point method ...")
        feasible_point_sets, fixed_reaction_ind_all, indp_spec_ind_dec, decision_vector_bounds = self.__feasible_point_method(bounds, iterations, seed, print_flag)
        print("Feasible point method has finished.")
        bounds_2 = decision_vector_bounds

        full_set_of_values = numpy.zeros(len(bounds), dtype=numpy.float64)

        inputs = numpy.zeros(len(self.__vars_for_lam_fixed_reactions), dtype=numpy.float64)

        input_values = numpy.zeros(len(self.__lambda_variables), dtype=numpy.float64)

        det_point_sets = []
        det_point_sets_fun = []
        smallest_value = numpy.float(1e8)
        print("Starting determinant optimization ...")
        for i in range(len(feasible_point_sets)):

            with numpy.errstate(divide='ignore', invalid='ignore'):

                result = scipy.optimize.dual_annealing(self.__obj_func, bounds=bounds_2, args=(self.__lambda_det_jac, bounds,
                                                                                full_set_of_values, fixed_reaction_ind_all,
                                                                                self.__cons_laws_lamb, indp_spec_ind_dec,
                                                                                inputs, input_values), x0=feasible_point_sets[i], seed=seed,
                                                       local_search_options={'method': "Nelder-Mead"}, maxiter=dual_annealing_iters)

                if print_flag:
                    print("Global function value: " + str(result.fun))
                    print("Decision vector used: ")
                    print(result.x)

                if abs(result.fun) > numpy.float64(1e-100):

                    result1 = scipy.optimize.minimize(self.__obj_func, result.x, args=(self.__lambda_det_jac, bounds,
                                                                                       full_set_of_values, fixed_reaction_ind_all,
                                                                                       self.__cons_laws_lamb, indp_spec_ind_dec,
                                                                                       inputs, input_values),
                                                     method='Nelder-Mead', tol=1e-16)

                    if print_flag:
                        print("Local function value: " + str(result1.fun))
                        print("Decision vector used: ")
                        print(result1.x)

                    if smallest_value > result1.fun:
                        smallest_value = result1.fun

                    if abs(result1.fun) <= numpy.finfo(float).eps:
                        det_point_sets.append(result1.x)
                        det_point_sets_fun.append(result1.fun)

                else:
                    if smallest_value > result.fun:
                        smallest_value = result.fun
                    det_point_sets.append(result.x)
                    det_point_sets_fun.append(result.fun)

                if print_flag:
                    print("")

        print("Determinant optimization has finished.")

        params = self.__final_check(det_point_sets, bounds, full_set_of_values, fixed_reaction_ind_all,
                                    self.__cons_laws_lamb, indp_spec_ind_dec, inputs, input_values)

        end_t = time.time()
        elapsed = end_t - start_t
        print("Determinant optimization has finished.")
        print("Elapsed time for optimization: " + str(elapsed))
        print("\nSmallest value achieved by objective function: " + str(smallest_value) + "\n")

        return params, det_point_sets_fun

    def __final_check(self, det_point_sets, boundz, full_set_of_values, reaction_ind,
                      cons_laws_lamb, indp_spec_ind, inputs, input_values):

        params = []
        for x in det_point_sets:

            # filling in values for non-fixed reactions
            inputs[0:self.__R - len(reaction_ind)] = x[0:self.__R - len(reaction_ind)]

            # filling in values for independent species
            count = self.__R - len(reaction_ind)
            for i in indp_spec_ind:
                inputs[count] = x[i]
                count += 1

            # filling in values for conservation laws
            for i in cons_laws_lamb:
                inputs[count] = i(*tuple(x[self.__R - len(reaction_ind):]))
                count += 1

            # filling in all species values
            full_set_of_values[self.__R:] = x[self.__R - len(reaction_ind):]

            # filling in non-fixed reaction values
            count = 0
            for i in range(self.__R):
                if i not in reaction_ind:
                    full_set_of_values[i] = x[count]
                    count += 1

            # filling in fixed reaction values
            full_set_of_values[reaction_ind] = numpy.array([i(*tuple(inputs)) for i in self.__lambda_fixed_reactions])

            # filing in reactions
            input_values[0: self.__R] = full_set_of_values[0:self.__R]

            # filling in independent species
            input_values[self.__R: self.__R + len(indp_spec_ind)] = x[indp_spec_ind]

            # filling in conservation laws
            input_values[self.__R + len(indp_spec_ind):] = inputs[self.__R - len(reaction_ind) + len(indp_spec_ind):]

            check = self.__saddle_node_bifurcation(input_values)
            # print("saddle node check:")
            # print(check)
            # print("")
            if True: #check:
                params.append(input_values.flatten())

        return params


    def __saddle_node_bifurcation(self, result_x):

        its_zero = numpy.float(1e-14) # necessary to compare against zero

        jac_w_vals = self.__lambda_jec(*tuple(result_x)).astype(numpy.float64)
        jac_nullspace = scipy.linalg.null_space(jac_w_vals).astype(numpy.float64)

        # check for one dimensional nullspace
        if jac_nullspace.shape[1] == 1:

            k_vec = jac_nullspace.astype(numpy.float64)

            # check for nonzero nullspace
            if any([abs(i) > its_zero for i in k_vec]):

                eigs = numpy.linalg.eigvals(jac_w_vals).real.astype(numpy.float64)
                ind_zero_sum = sum([int(1) for i in eigs if abs(i) <= its_zero])

                # check for only one zero eigenvalue
                if ind_zero_sum == 1:

                    # check for steady state
                    if all([i(*tuple(result_x)) <= its_zero for i in self.__lambda_indp_odes]):

                        jac_eps_val = self.__lambda_jac_eps(*tuple(result_x)).astype(numpy.float64)

                        k_vec_vals = [numpy.float64(i) for i in k_vec]
                        D2_g_u_w_val = self.__lambda_D2_g_u_w(*tuple(list(result_x) + k_vec_vals)).astype(numpy.float64)

                        Duu_kk = numpy.dot(D2_g_u_w_val, k_vec).astype(numpy.float64)

                        # check for nonzero f_eps and f_uu
                        if any([abs(i) > its_zero for i in jac_eps_val]) \
                                and any([abs(i) > its_zero for i in Duu_kk]):

                            aug_eps = numpy.column_stack((jac_w_vals, jac_eps_val))
                            aug_Duu = numpy.column_stack((jac_w_vals, Duu_kk))

                            jac_w_vals_rank = numpy.linalg.matrix_rank(jac_w_vals)
                            aug_eps_rank = numpy.linalg.matrix_rank(aug_eps)
                            aug_Duu_rank = numpy.linalg.matrix_rank(aug_Duu)

                            # check that f_eps and f_uu are not in the range of f_u
                            if (jac_w_vals_rank < aug_eps_rank) and (jac_w_vals_rank < aug_Duu_rank):
                                return True
                            else:
                                return False
                        else:
                            return False
                    else:
                        return False
                else:
                    return False
            else:
                return False
        else:
            return False

    def run_greedy_continuity_analysis(self, species=None, parameters=None, dir_path="./num_cont_graphs",
                                       print_lbls_flag=False, auto_parameters=None):

        # setting default values for AUTO
        if 'NMX' not in auto_parameters.keys():
            auto_parameters['NMX'] = 10000

        if 'ITMX' not in auto_parameters.keys():
            auto_parameters['ITMX'] = 100

        # making the directory if it doesn't exist
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)

        species_num = [str(i) for i in self.__indp_species].index(species) + 1

        species_y = str(self.__indp_species[species_num-1])

        multistable_param_ind, important_info, plot_specifications = BistabilityFinder.run_greedy_continuity_analysis \
            (species_num, parameters, self.__initialize_ant_string, self.__finalize_ant_string, species_y, dir_path,
             print_lbls_flag, auto_parameters)

        self.__important_info += important_info

        return multistable_param_ind, plot_specifications

    def __initialize_ant_string(self, species_num, pcp_x):

        indp_odes_str = []
        indp_odes_str.append(str(self.__indp_system_subs[species_num-1]))
        [indp_odes_str.append(str(self.__indp_system_subs[i])) for i in range(len(self.__indp_species)) if i != species_num-1]

        for i in range(len(indp_odes_str)):
            indp_odes_str[i] = indp_odes_str[i].replace('**', '^')

        ind_spec_conc = []
        ind_spec_conc.append(self.__indp_species[species_num-1])
        [ind_spec_conc.append(str(i)) for i in self.__indp_species if i != self.__indp_species[species_num-1]]

        # building the string of ODEs in Antimony syntax
        ode_str = ''
        for i in range(len(ind_spec_conc)):
            ode_str += 'J' + str(i) + ': -> ' + str(ind_spec_conc[i]) + '; ' + indp_odes_str[i] + ';'

        return ode_str, pcp_x

    def __finalize_ant_string(self, result_x, ode_str):

        vars_to_initialize = [str(i) for i in self.__lambda_variables]
        var_vals = result_x

        ant_str = ode_str
        for i in range(len(vars_to_initialize)):
            ant_str += str(vars_to_initialize[i]) + ' = ' + str(var_vals[i]) + ';'

        return ant_str

    def get_optimization_bounds(self):

        graph_edges = self.__cgraph.get_g_edges()
        dec_vec_var_def = []
        for i in self.__sympy_reactions + self.__sympy_species:

            if i in self.__sympy_species:
                dec_vec_var_def.append("concentration")
            elif i in self.__sympy_reactions:

                ind = self.__sympy_reactions.index(i)
                reaction = graph_edges[ind]

                reaction_type = self.__cgraph.get_graph().edges[reaction]['type']
                dec_vec_var_def.append(reaction_type)

                if reaction_type is None:
                    output_statement = "The reaction type of reaction " + self.__cgraph.get_graph().edges[reaction][
                        'label'] \
                                       + " could not be identified as it does not fit any biological criteria " + \
                                       "established. \n" + "You must enter bounds manually for this reaction! \n"
                    print(output_statement)

        bounds = [self.get_physiological_range(i) for i in dec_vec_var_def]

        return bounds

    @staticmethod
    def get_physiological_range(for_what=None):
        """
        Molecule concentrations: 5e-13 ... 5e-7 M (basically from 1 per cell to 1e6 per cell)
        Complex formation (bimolecular association): 1e4 ... 1e8  M^-1s^-1
        Example A + B -> [AB]  (or first step in enzymatic catalysis mechanism)
        Complex dissociation: 1e-5 ... 1e-3 s^-1
        Example [AB] -> A + B  (e.g enzyme substrat dissociation)
        Catalysis (equivalent to kcat in Michaelis-Menten): 1e-3 ... 1e0 s^-1
        Examples: [AB] -> A + C  or  X -> Y  (e.g. final step in enzymatic catalysis)
        """
        valid = {"concentration", "complex formation", "complex dissociation", "catalysis"}
        if for_what is None:
            warnings.warn("Please provide the argument what the range should be provided for.")
            return None
        if for_what not in valid:
            raise ValueError("for_what argument must be one of %r." % valid)
        if for_what == "concentration":
            return 5e-13, 5e-7
        if for_what == "complex formation":
            return 1e4, 1e8
        if for_what == "complex dissociation":
            return 1e-5, 1e-3
        if for_what == "catalysis":
            return 1e-3, 1e0
    # @staticmethod
    # def get_physiological_range(for_what=None):
    #     """Obtains physiological ranges.
    #
    #     Parameters
    #     -----------
    #     for_what: string
    #         Accepted values: "concentration", "complex formation", "complex dissociation", "catalysis", or "flux"
    #
    #     Returns
    #     --------
    #     concentration: tuple
    #         (5e-1,5e5) pM
    #     complex formation: tuple
    #         (1e-8,1e-4)  pM^-1s^-1
    #     complex dissociation: tuple
    #         (1e-5,1e-3) s^-1
    #     catalysis: tuple
    #         (1e-3,1) s^-1
    #     flux: tuple
    #         (0, 55) M s^-1
    #
    #
    #     Example
    #     --------
    #     >>> import crnt4sbml
    #     >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
    #     >>> network.get_physiological_range("concentration")
    #     """
    #     valid = {"concentration", "complex formation", "complex dissociation", "catalysis", "flux"}
    #     if for_what is None:
    #         warnings.warn("Please provide the argument what the range should be provided for.")
    #         return None
    #     if for_what not in valid:
    #         raise ValueError("for_what argument must be one of %r." % valid)
    #     if for_what == "concentration":
    #         return (5e-1, 5e5)
    #     if for_what == "complex formation":
    #         return (1e-8, 1e-4)
    #     if for_what == "complex dissociation":
    #         return (1e-5, 1e-3)
    #     if for_what == "catalysis":
    #         return (1e-3, 1e0)
    #     if for_what == "flux":
    #         return (0, 55)  # e12)

    def get_input_vector(self):
        return self.__sympy_species + self.__sympy_reactions