# general appraoch for any network that conserves mass 
# Subclass of BistabilityFinder and BistabilityAnalysis
import sympy
import sys
import itertools
import numpy
import scipy.optimize
import scipy.linalg
import warnings
import time
import math
import os
from .bistability_finder import BistabilityFinder
from .bistability_analysis import BistabilityAnalysis


class GeneralApproach(BistabilityFinder, BistabilityAnalysis):
    """
    Class for constructing a more general approach to bistability detection for systems with mass action kinetics.
    """
    def __init__(self, cgraph, get_physiological_range):
        """
        Initialization of GeneralApproach class.

        See also
        ---------
        crnt4sbml.CRNT.get_general_approach
        """
        self.__cgraph = cgraph
        self.__get_physiological_range = get_physiological_range

        if not self.__cgraph.get_dim_equilibrium_manifold() > 0:
            print("# of species - rank(S) is not greater than zero.")
            print("This implies that there are no mass conservation laws.")
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

        self.__create_BT_matrix()

        self.__cons_laws_sympy = None
        self.__indp_species = None
        self.__modified_ode_system = None
        self.__species_mod_system = None
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
        self.__comm = None
        self.__my_rank = None
        self.__num_cores = None
        self.__method = "GeneralApproach"

        self.__create_conservation_laws()

    def initialize_general_approach(self, signal=None, response=None, fix_reactions=False):
        """
        Function for initializing the necessary variables for the general approach.

        Parameters
        -----------
            signal: String
                A string stating the conservation law that is the x-axis of the bifurcation diagram.
            response: String
                A string stating the species that is the y-axis of the bifurcation diagram.
            fix_reactions: bool
                A bool that determines if a steady state is enforced by fixing the reactions. See :ref:`gen-app-label`
                for specific details.

        Examples
        ---------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> GA = network.get_general_approach()
        >>> signal = "C1"
        >>> response = "s1"
        >>> GA.initialize_general_approach(signal=signal, response=response)

        See also :ref:`quickstart-gen-app-label` and :ref:`gen-app-label`.
        """
        if signal in ['C' + str(i + 1) for i in range(len(self.__cons_laws_sympy))] \
                and response in self.__cgraph.get_species():

            self.__signal = signal
            self.__response = response
            self.__fix_reactions = fix_reactions

            self.__create_modified_ode_system()
            self.__create_variables_for_optimization()

        else:
            print("The provided signal, response, or both are not properly defined.")
            print("This may be the result of the values not being provided or the signal is not a conservation law or")
            print("the response is not a species concentration.")

    def __create_BT_matrix(self):
        # (bt = NullSpace[Transpose[S]]) // MatrixForm
        the_null_space = self.__S_mat.T.nullspace()

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

        self.__cons_laws_sympy_lamb = [sympy.utilities.lambdify(self.__sympy_species, self.__cons_laws_sympy[i]) for i in range(len(self.__cons_laws_sympy))]

    def __is_list_empty(self, inlist):
        if isinstance(inlist, list):  # Is a list
            return all(map(self.__is_list_empty, inlist))
        return False  # Not a list

    def __create_modified_ode_system(self):

        conservation_constants = ['C' + str(i + 1) for i in range(len(self.__cons_laws_sympy))]

        self.__signal_index = conservation_constants.index(self.__signal)

        cons_laws_sympy_eq = [sympy.Eq(sympy.Symbol('C' + str(i + 1), real=True), self.__cons_laws_sympy[i])
                              for i in range(len(self.__cons_laws_sympy))]

        elements = self.__cons_laws_sympy[self.__signal_index].atoms()

        possible_species = [i for i in self.__sympy_species if i in elements]

        dep_species = [i for i in possible_species if i != self.__response][0]

        dep_species_expression = sympy.solve(cons_laws_sympy_eq[self.__signal_index], dep_species)[0]

        dep_species_index = self.__sympy_species.index(dep_species)

        self.__modified_ode_system = sympy.Matrix([self.__full_system[i] for i in range(self.__N) if i != dep_species_index])

        # substitute in expression for dependent species
        for i in range(self.__N-1):
            self.__modified_ode_system[i] = self.__modified_ode_system[i].subs(dep_species, dep_species_expression)

        self.__species_mod_system = [i for i in self.__sympy_species if i != dep_species]

        self.__construct_independent_system()
        self.__construct_independent_system_with_substitution()

        self.__x_bar = self.__sympy_reactions + self.__sympy_species
        self.__lagrangian_vars = self.__x_bar + [sympy.Symbol('C' + str(i + 1), real=True) for i in
                                                 range(len(self.__cons_laws_sympy))]

        if self.__fix_reactions:
            self.__enforce_steady_state()
            self.__fixed_reaction_indices = [i for i in range(len(self.__lagrangian_vars)) if
                                             self.__lagrangian_vars[i] in self.__fixed_reactions]

    @staticmethod
    def __unique_everseen(iterable, key=None):
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
        possible_fixed_reactions = list(self.__unique_everseen(temp, key=frozenset))

        return possible_fixed_reactions

    def __enforce_steady_state(self):

        self.__fixed_reactions = []
        a, vals = self.__S_mat.rref()

        self.__fixed_reactions = [self.__sympy_reactions[i] for i in vals]
        a, b = sympy.linear_eq_to_matrix(self.__indp_system, self.__fixed_reactions)
        temp_solution_tuple = list(sympy.linsolve((a, b), self.__fixed_reactions))[0]

        no_zero_values_flag = True

        # looking for other fixed reaction sets if no solution was found
        if sympy.S.Zero in temp_solution_tuple:
            print("")
            print("At least one fixed reaction found to be zero. \nAttempting to find nonzero values for all fixed reactions.")

            no_zero_values_flag = False
            reactions_found_to_be_zero = []

            zero_indices = [i for i in range(len(temp_solution_tuple)) if temp_solution_tuple[i] == sympy.S.Zero]
            [reactions_found_to_be_zero.append(self.__fixed_reactions[i]) for i in zero_indices]

            comb = self.__construct_possible_fixed_reactions()
            zeros = sympy.zeros(len(self.__indp_system), 1)
            for i in comb:
                i.sort(key=lambda j: self.__sympy_reactions.index(j))
                a, b = sympy.linear_eq_to_matrix(self.__indp_system, i)
                rank_a = a.rank()
                aug = a.row_join(zeros)
                rank_aug = aug.rank()
                self.__fixed_reactions = i
                if rank_a == rank_aug and rank_a == a.shape[1]:
                    self.__fixed_reactions = i
                    temp_solution_tuple2 = list(sympy.linsolve((a, b), self.__fixed_reactions))[0]
                    if sympy.S.Zero not in temp_solution_tuple2:
                        no_zero_values_flag = True
                        break
                    else:
                        zero_indices = [ii for ii in range(len(temp_solution_tuple2)) if
                                        temp_solution_tuple2[ii] == sympy.S.Zero]
                        [reactions_found_to_be_zero.append(i[ii]) for ii in zero_indices if i[ii] not in reactions_found_to_be_zero]
            print("Done searching for nonzero values for all fixed reactions.")
            print("")

        if no_zero_values_flag:
            self.__soln_to_fixed_reactions = list(temp_solution_tuple)
        else:
            print("There was no solution which resulted in all fixed reaction values being nonzero.")
            print(f"Consider removing one or all of the following reactions {reactions_found_to_be_zero}. \n")
            self.__soln_to_fixed_reactions = list(temp_solution_tuple)

        for i in range(len(self.__soln_to_fixed_reactions)):
            for j in range(len(self.__replacements)):
                self.__soln_to_fixed_reactions[i] = self.__soln_to_fixed_reactions[i].subs(self.__replacements[j][0], self.__replacements[j][1])

        self.__vars_for_lam_fixed_reactions = [i for i in self.__lagrangian_vars if i not in self.__fixed_reactions]
        self.__lambda_fixed_reactions = [sympy.utilities.lambdify(self.__vars_for_lam_fixed_reactions, i) for i in
                                         self.__soln_to_fixed_reactions]

    def __create_variables_for_optimization(self):

        self.__jac_mod_system = self.__indp_system_subs.jacobian(sympy.Matrix(self.__indp_species))

        self.__det_jac = self.__jac_mod_system.det(method='lu')

        self.__build_lagrangian()

        self.__build_objective_function()

    def __build_lagrangian(self):

        if self.__fix_reactions:
            self.__lagrangian = self.__det_jac**2
        else:
            second_term = sympy.S.Zero
            for i in self.__indp_system_subs:
                second_term += i**2

            self.__lagrangian = self.__det_jac**2 + second_term

            # enforcing a steady state
            self.__steady_state_func = second_term

            self.__steady_state_lambda = sympy.utilities.lambdify(self.__lagrangian_vars, self.__steady_state_func)

        self.__lambda_lagrangian = sympy.utilities.lambdify(self.__lagrangian_vars, self.__lagrangian)

    def __build_objective_function(self):

        self.__obj_func_h = sympy.S.Zero

        self.__obj_func_h = self.__lagrangian

        self.__lambda_obj_func_h = self.__lambda_lagrangian

    def __obj_func_fixed(self, x):

        count0 = 0
        for i in range(len(self.__x_full) - len(self.__cons_laws_sympy_lamb)):
            if i not in self.__fixed_reaction_indices:
                self.__x_full[i] = x[count0]
                count0 += 1

        count = 0
        for i in self.__cons_laws_sympy_lamb:
            self.__x_full[len(self.__x_bar) + count] = i(*tuple(x[self.__R-len(self.__fixed_reaction_indices):
                                                                  self.__R-len(self.__fixed_reaction_indices) + self.__N]))
            count += 1

        for i in range(len(self.__lambda_fixed_reactions)):
            self.__x_full[self.__fixed_reaction_indices[i]] = self.__lambda_fixed_reactions[i](*tuple([self.__x_full[i] for i in range(len(self.__x_full)) if i not in self.__fixed_reaction_indices]))

        fun = self.__lambda_obj_func_h(*tuple(self.__x_full))

        # bounds of the optimization problem
        sumval = numpy.float64(0.0)
        for i in range(len(self.__x_full) - len(self.__cons_laws_sympy_lamb)):
            sumval += numpy.maximum(numpy.float64(0.0), numpy.float64(self.__bounds[i][0]) - self.__x_full[i])
            sumval += numpy.maximum(numpy.float64(0.0), self.__x_full[i] - numpy.float64(self.__bounds[i][1]))

        fun += sumval

        # constraints of the optimization problem
        if self.__full_constraints:
            for i in range(len(self.__full_constraints)):
                fun += self.__full_constraints[i][0](self.__x_full[0:len(self.__x_full) - len(self.__cons_laws_sympy_lamb)], self.__full_constraints[i][1])

        if numpy.isnan(fun):
            return numpy.Inf
        else:
            return fun

    def __obj_func(self, x):

        self.__x_full[0:len(self.__x_bar)] = x

        count = 0
        for i in self.__cons_laws_sympy_lamb:
            self.__x_full[len(self.__x_bar) + count] = i(*tuple(x[self.__R:self.__R + self.__N]))
            count += 1

        fun = self.__lambda_obj_func_h(*tuple(self.__x_full))

        # bounds of the optimization problem
        sumval = numpy.float64(0.0)
        for i in range(len(x)):
            sumval += numpy.maximum(numpy.float64(0.0), numpy.float64(self.__bounds[i][0]) - x[i])
            sumval += numpy.maximum(numpy.float64(0.0), x[i] - numpy.float64(self.__bounds[i][1]))

        fun += sumval

        # constraints of the optimization problem
        if self.__full_constraints:
            for i in range(len(self.__full_constraints)):
                fun += self.__full_constraints[i][0](self.__x_full[0:len(self.__x_full) - len(self.__cons_laws_sympy_lamb)], self.__full_constraints[i][1])

        if numpy.isnan(fun):
            return numpy.Inf
        else:
            return fun

    def run_optimization(self, bounds=None, iterations=10, seed=0, print_flag=False, dual_annealing_iters=1000,
                         confidence_level_flag=False, change_in_rel_error=1e-2, constraints=None, parallel_flag=False):
        """
        Function for running the optimization problem for the general approach.

        Parameters
        -----------
            bounds: list of tuples
                A list defining the lower and upper bounds for each variable in the input vector. See
                :func:`crnt4sbml.GeneralApproach.get_input_vector`.
            iterations: int
                The number of iterations to run the multistart method.
            seed: int
                Seed for the random number generator. None should be used if a random generation is desired.
            print_flag: bool
                Should be set to True if the user wants the objective function values found in the optimization problem
                and False otherwise.
            dual_annealing_iters: integer
                The number of iterations that should be ran for dual annealing routine in optimization.
            confidence_level_flag: bool
                If True a confidence level for the objective function will be given.
            change_in_rel_error: float
                The maximum relative error that should be allowed to consider :math:`f_k` in the neighborhood
                of :math:`\widetilde{f}`.
            constraints: list of dictionaries
                Each dictionary is of the form {'type': '...', 'fun': lambda x:  ... }, where 'type' can be set to
                'ineq' or 'eq' and the lambda function to be defined by the user. The 'ineq' refers to an
                inequality constraint c(x) with c(x) <= 0. For the lambda function the input x refers to the input vector
                of the optimization routine. See :ref:`gen-app-label` for further details.
            parallel_flag: bool
                If set to True a parallel version of the optimization routine is ran. If False, a serial version of the
                optimization routine is ran. See :ref:`parallel-gen-app-label`.
        Returns
        --------
        params_for_global_min: list of numpy arrays
            A list of numpy arrays that correspond to the input vectors of the problem.
        obj_fun_val_for_params: list of floats
            A list of objective function values produced by the corresponding input vectors in params_for_global_min.

        Examples
        ---------
        See :ref:`quickstart-gen-app-label` and :ref:`gen-app-label`.
        """

        self.__initialize_optimization_variables(bounds, iterations, seed, print_flag, dual_annealing_iters,
                                                 confidence_level_flag, change_in_rel_error, constraints, parallel_flag)

        params_for_global_min, obj_fun_val_for_params, self.__important_info = self._BistabilityFinder__parent_run_optimization()

        self.__my_rank = self._BistabilityFinder__my_rank
        self.__comm = self._BistabilityFinder__comm

        return params_for_global_min, obj_fun_val_for_params

    def __initialize_optimization_variables(self, bounds, iterations, seed, print_flag, dual_annealing_iters,
                                            confidence_level_flag, change_in_rel_error, constraints, parallel_flag):

        self.__bounds = bounds
        self.__iterations = iterations
        self.__seed = seed
        self.__print_flag = print_flag
        self.__dual_annealing_iters = dual_annealing_iters
        self.__confidence_level_flag = confidence_level_flag
        self.__change_in_rel_error = change_in_rel_error
        self.__constraints = constraints
        self.__parallel_flag = parallel_flag
        self.__x_full = None
        self.__full_constraints = None
        self.__temp_bounds = None

    def __run_global_optimization_routine(self, initial_x):

        if self.__fix_reactions:
            result = scipy.optimize.dual_annealing(self.__obj_func_fixed, bounds=self.__temp_bounds, x0=initial_x,
                                                   seed=self.__seed, local_search_options={'method': "Nelder-Mead"},
                                                   maxiter=self.__dual_annealing_iters)
        else:
            result = scipy.optimize.dual_annealing(self.__obj_func, bounds=self.__bounds, x0=initial_x,
                                                   seed=self.__seed, local_search_options={'method': "Nelder-Mead"},
                                                   maxiter=self.__dual_annealing_iters)
        return result

    def __run_local_optimization_routine(self, initial_x):

        if self.__fix_reactions:
            result1 = scipy.optimize.minimize(self.__obj_func_fixed, initial_x, method='Nelder-Mead', tol=1e-16)
        else:
            result1 = scipy.optimize.minimize(self.__obj_func, initial_x, method='Nelder-Mead', tol=1e-16)

        return result1

    def __create_final_points(self, x):

        if self.__fix_reactions:

            count0 = 0
            for i in range(len(self.__x_full) - len(self.__cons_laws_sympy_lamb)):
                if i not in self.__fixed_reaction_indices:
                    self.__x_full[i] = x[count0]
                    count0 += 1

            count = 0
            for i in self.__cons_laws_sympy_lamb:
                self.__x_full[len(self.__x_bar) + count] = i(*tuple(x[self.__R - len(
                    self.__fixed_reaction_indices):self.__R - len(self.__fixed_reaction_indices) + self.__N]))
                count += 1

            for i in range(len(self.__lambda_fixed_reactions)):
                self.__x_full[self.__fixed_reaction_indices[i]] = self.__lambda_fixed_reactions[i](
                    *tuple([self.__x_full[i] for i in range(len(self.__x_full)) if i not in self.__fixed_reaction_indices]))

            return numpy.array(list(self.__x_full[0:len(self.__x_bar)]))

        else:
            return numpy.array(list(x[:]))

    def run_greedy_continuity_analysis(self, species=None, parameters=None, dir_path="./num_cont_graphs",
                                       print_lbls_flag=False, auto_parameters=None, plot_labels=None):

        """
        Function for running the greedy numerical continuation and bistability analysis portions of the general
        approach. This routine uses the initial value of the principal continuation parameter to construct AUTO
        parameters and then tests varying fixed step sizes for the continuation problem. Note that this routine may
        produce jagged or missing sections in the plots provided. To produce better plots one should use the information
        provided by this routine to run :func:`crnt4sbml.GeneralApproach.run_continuity_analysis`.

        Note: A parallel version of this routine is not available.

        Parameters
        ------------
            species: string
                A string stating the species that is the y-axis of the bifurcation diagram.
            parameters: list of numpy arrays
                A list of numpy arrays corresponding to the decision vectors that produce a small objective function
                value.
            dir_path: string
                A string stating the path where the bifurcation diagrams should be saved.
            print_lbls_flag: bool
                If True the routine will print the special points found by AUTO 2000 and False will not print any
                special points.
            auto_parameters: dict
                Dictionary defining the parameters for the AUTO 2000 run. Please note that only the
                PrincipalContinuationParameter in this dictionary should be defined, no other AUTO parameters should
                be set. For more information on these parameters refer to :download:`AUTO parameters <../auto2000_input.pdf>`.
            plot_labels: list of strings
                A list of strings defining the labels for the x-axis, y-axis, and title. Where the first element
                is the label for x-axis, second is the y-axis label, and the last element is the title label. If
                you would like to use the default settings for some of the labels, simply provide None for that
                element.
        Returns
        ---------
            multistable_param_ind: list of integers
                A list of those indices in 'parameters' that produce multistable plots.
            plot_specifications: list of lists
                A list whose elements correspond to the plot specifications of each element in multistable_param_ind.
                Each element is a list where the first element specifies the range used for the x-axis, the second
                element is the range for the y-axis, and the last element provides the x-y values and special point label
                for each special point in the plot.

        Example
        ---------
        See :ref:`quickstart-gen-app-label` and :ref:`gen-app-label`.
        """
        if self.__comm is not None:
            if self.__my_rank == 0:
                self.__initialize_continuity_analysis(species, parameters, dir_path, print_lbls_flag, auto_parameters,
                                                      plot_labels)
                multistable_param_ind, important_info, plot_specifications = self._BistabilityAnalysis__parent_run_greedy_continuity_analysis()
            else:
                important_info = ''
                multistable_param_ind = []
                plot_specifications = []
            self.__comm.Barrier()
        else:
            self.__initialize_continuity_analysis(species, parameters, dir_path, print_lbls_flag, auto_parameters,
                                                  plot_labels)
            multistable_param_ind, important_info, plot_specifications = self._BistabilityAnalysis__parent_run_greedy_continuity_analysis()

        self.__important_info += important_info

        return multistable_param_ind, plot_specifications

    def run_continuity_analysis(self, species=None, parameters=None, dir_path="./num_cont_graphs",
                                       print_lbls_flag=False, auto_parameters=None, plot_labels=None):

        """
        Function for running the numerical continuation and bistability analysis portions of the general
        approach.

        Note: A parallel version of this routine is not available.

        Parameters
        ------------
            species: string
                A string stating the species that is the y-axis of the bifurcation diagram.
            parameters: list of numpy arrays
                A list of numpy arrays corresponding to the decision vectors that produce a small objective function
                value.
            dir_path: string
                A string stating the path where the bifurcation diagrams should be saved.
            print_lbls_flag: bool
                If True the routine will print the special points found by AUTO 2000 and False will not print any
                special points.
            auto_parameters: dict
                Dictionary defining the parameters for the AUTO 2000 run. Please note that one should **not** set
                'SBML' or 'ScanDirection' in these parameters as these are automatically assigned. It is absolutely
                necessary to set PrincipalContinuationParameter in this dictionary. For more information on these
                parameters refer to :download:`AUTO parameters <../auto2000_input.pdf>`. 'NMX' will default to
                10000 and 'ITMX' to 100.
            plot_labels: list of strings
                A list of strings defining the labels for the x-axis, y-axis, and title. Where the first element
                is the label for x-axis, second is the y-axis label, and the last element is the title label. If
                you would like to use the default settings for some of the labels, simply provide None for that
                element.
        Returns
        ---------
            multistable_param_ind: list of integers
                A list of those indices in 'parameters' that produce multistable plots.
            plot_specifications: list of lists
                A list whose elements correspond to the plot specifications of each element in multistable_param_ind.
                Each element is a list where the first element specifies the range used for the x-axis, the second
                element is the range for the y-axis, and the last element provides the x-y values and special point label
                for each special point in the plot.

        Example
        ---------
        See :ref:`quickstart-gen-app-label` and :ref:`gen-app-label`..
        """
        if self.__comm is not None:
            if self.__my_rank == 0:
                self.__initialize_continuity_analysis(species, parameters, dir_path, print_lbls_flag, auto_parameters,
                                                      plot_labels)
                multistable_param_ind, important_info, plot_specifications = self._BistabilityAnalysis__parent_run_continuity_analysis()
            else:
                important_info = ''
                multistable_param_ind = []
                plot_specifications = []
            self.__comm.Barrier()
        else:
            self.__initialize_continuity_analysis(species, parameters, dir_path, print_lbls_flag, auto_parameters,
                                                  plot_labels)
            multistable_param_ind, important_info, plot_specifications = self._BistabilityAnalysis__parent_run_continuity_analysis()

        self.__important_info += important_info

        return multistable_param_ind, plot_specifications

    def __initialize_continuity_analysis(self, species, parameters, dir_path, print_lbls_flag, auto_parameters, plot_labels):

        self.__parameters = parameters
        self.__dir_path = dir_path
        self.__print_lbls_flag = print_lbls_flag
        self.__auto_parameters = auto_parameters
        self.__plot_labels = plot_labels

        if self.__comm is not None:

            print("")
            print("A parallel version of numerical continuation is not available.")
            print("Numerical continuation will be ran using only one core.")
            print("For your convenience, the provided parameters have been saved in the current directory under the name params.npy.")
            numpy.save('./params.npy', parameters)

        # setting default values for AUTO
        if 'NMX' not in self.__auto_parameters.keys():
            self.__auto_parameters['NMX'] = 10000

        if 'ITMX' not in self.__auto_parameters.keys():
            self.__auto_parameters['ITMX'] = 100

        # making the directory if it doesn't exist
        if not os.path.isdir(self.__dir_path):
            os.mkdir(self.__dir_path)

        self.__species_num = [str(i) for i in self.__indp_species].index(species) + 1

        self.__species_y = str(self.__indp_species[self.__species_num - 1])


    def run_direct_simulation(self, params_for_global_min=None, dir_path="./dir_sim_graphs", change_in_relative_error=1e-6,
                              parallel_flag=False, print_flag=False, left_multiplier=0.5, right_multiplier=0.5):

        """
        Function for running direct simulation to conduct bistability analysis of the general approach.

        Note: This routine is more expensive than the numerical continuation routines, but can provide solutions
        when the Jacobian of the ODE system is always singular. A parallel version of this routine is available.

        Parameters
        ------------
            params_for_global_min: list of numpy arrays
                A list of numpy arrays corresponding to the input vectors that produce a small objective function
                value.
            dir_path: string
                A string stating the path where the bifurcation diagrams should be saved.
            change_in_relative_error: float
                A float value that determines how small the relative error should be in order for the solution of the
                ODE system to be considered at a steady state. Note: a smaller value will run faster, but may produce
                an ODE system that is not at a steady state.
            parallel_flag: bool
                If set to True a parallel version of direct simulation is ran. If False, a serial version of the
                routine is ran. See :ref:`parallel-gen-app-label` for further information.
            print_flag: bool
                If set to True information about the direct simulation routine will be printed. If False, no output
                will be provided.
            left_multiplier: float
                A float value that determines the percentage of the signal that will be searched to the left of the signal
                value. For example, the lowerbound for the signal range will be signal_value - signal_value*left_multiplier.
            right_multiplier: float
                A float value that determines the percentage of the signal that will be searched to the right of the signal
                value. For example, the upperbound for the signal range will be signal_value + signal_value*right_multiplier.

        Example
        ---------
        See :ref:`gen-app-label`.
        """
        self.__initialize_direct_simulation(params_for_global_min, dir_path, change_in_relative_error, parallel_flag,
                                            print_flag, left_multiplier, right_multiplier)

        self._BistabilityAnalysis__parent_run_direct_simulation()

        self.__my_rank = self._BistabilityAnalysis__my_rank
        self.__comm = self._BistabilityAnalysis__comm

    def __initialize_direct_simulation(self, params_for_global_min, dir_path, change_in_relative_error, parallel_flag,
                                       print_flag, left_multiplier, right_multiplier):

        self.__parameters = params_for_global_min
        self.__dir_path = dir_path
        self.__change_in_relative_error = change_in_relative_error
        self.__parallel_flag = parallel_flag
        self.__dir_sim_print_flag = print_flag
        self.__left_multiplier = left_multiplier
        self.__right_multiplier = right_multiplier
        self.__concentration_funs = None

        lambda_inputs = self.__sympy_reactions + self.__sympy_species
        self.__ode_lambda_functions = [sympy.utilities.lambdify(lambda_inputs, self.__full_system[i]) for i in
                                       range(len(self.__full_system))]

        self.__jac_lambda_function = sympy.utilities.lambdify(lambda_inputs, self.__full_system.jacobian(self.__sympy_species))

    def __construct_independent_system(self):

        # forming ya matrix
        #ya = self.__Y_mat * self.__A_mat   # old approach
        ya, __ = sympy.linear_eq_to_matrix(self.__Y_mat * self.__A_mat*self.__psi_vec, self.__sympy_reactions) # new

        # finding how many rows are indep in ya
        __, vals = ya.T.rref()

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
        indp_odes_temp = ya[chosen_indp_indices, :] * sympy.Matrix(self.__sympy_reactions) # new approach

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

        vars_to_initialize = [str(i) for i in self.__x_bar] + ['C' + str(i + 1)
                                                               for i in range(len(self.__cons_laws_sympy))]

        conservation_vals = [self.__cons_laws_sympy_lamb[i](*tuple(result_x[self.__R:self.__R + self.__N]))
                             for i in range(len(self.__cons_laws_sympy))]

        var_vals = result_x

        ant_str = ode_str
        for i in range(len(self.__x_bar)):
            ant_str += str(vars_to_initialize[i]) + ' = ' + str(var_vals[i]) + ';'

        count = 0
        for i in range(len(self.__x_bar), len(vars_to_initialize)):
            ant_str += str(vars_to_initialize[i]) + ' = ' + str(conservation_vals[count]) + ';'
            count += 1

        if self.__print_lbls_flag:
            print(ant_str)

        return ant_str

    def get_optimization_bounds(self):
        """
        Returns a list of tuples that corresponds to the determined physiological bounds chosen for the problem. Each
        entry corresponds to the list provided by :func:`crnt4sbml.GeneralApproach.get_input_vector`.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> GA = network.get_general_approach()
        >>> signal = "C1"
        >>> response = "s1"
        >>> GA.initialize_general_approach(signal=signal, response=response)
        >>> GA.get_optimization_bounds()
        """
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

        bounds = [self.__get_physiological_range(i) for i in dec_vec_var_def]

        return bounds

    def generate_report(self):
        """
        Prints out helpful details constructed by :func:`crnt4sbml.GeneralApproach.run_optimization`,
        :func:`crnt4sbml.GeneralApproach.run_continuity_analysis`, and
        :func:`crnt4sbml.GeneralApproach.run_greedy_continuity_analysis`.

        Example
        --------
        See also :ref:`quickstart-gen-app-label` and :ref:`gen-app-label`
        """
        if self.__comm == None:
            print(self.__important_info)
        else:

            all_important_info = self.__comm.gather(self.__important_info, root=0)
            self.__comm.Barrier()

            if self.__my_rank == 0:
                for i in range(1, len(all_important_info)):
                    if all_important_info[i] != "":
                        print(all_important_info[i])
                print(self.__important_info)

    def get_input_vector(self):
        """
        Returns a list of SymPy variables that specifies the ordering of the reactions and species for which bounds
        need to be provided.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> GA = network.get_general_approach()
        >>> signal = "C1"
        >>> response = "s1"
        >>> GA.initialize_general_approach(signal=signal, response=response)
        >>> print(GA.get_input_vector())
        """
        return self.__sympy_reactions + self.__sympy_species

    def get_decision_vector(self):
        """
        Returns a list of SymPy variables that specifies the ordering of the reactions and species of the decision
        vector used in optimization. Note: this method should not be used to create bounds for the optimization
        routine, rather :func:`crnt4sbml.GeneralApproach.get_input_vector` should be used.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> GA = network.get_general_approach()
        >>> signal = "C1"
        >>> response = "s1"
        >>> GA.initialize_general_approach(signal=signal, response=response)
        >>> print(GA.get_decision_vector())
        """
        return self.__decision_vector

    def get_independent_odes_subs(self):
        """
        Returns a Sympy Matrix representing the independent ODE system with conservation laws substituted in. Each row
        corresponds to the ODE for the species corresponding to the list provided by
        :func:`crnt4sbml.GeneralApproach.get_independent_species`.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> GA = network.get_general_approach()
        >>> signal = "C1"
        >>> response = "s1"
        >>> GA.initialize_general_approach(signal=signal, response=response)
        >>> GA.get_independent_odes_subs()
        """
        return self.__indp_system_subs

    def get_independent_odes(self):
        """
        Returns a Sympy Matrix representing the independent ODE system without conservation laws substituted in. Each row
        corresponds to the ODE for the species corresponding to the list provided by :func:`crnt4sbml.GeneralApproach.get_independent_species`.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> GA = network.get_general_approach()
        >>> signal = "C1"
        >>> response = "s1"
        >>> GA.initialize_general_approach(signal=signal, response=response)
        >>> GA.get_independent_odes()
        """
        return self.__indp_system

    def get_independent_species(self):
        """
        Returns a list of SymPy variables that reflects the independent species chosen for the general approach.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> GA = network.get_general_approach()
        >>> signal = "C1"
        >>> response = "s1"
        >>> GA.initialize_general_approach(signal=signal, response=response)
        >>> GA.get_independent_species()
        """
        return self.__indp_species

    def get_fixed_reactions(self):
        """
        Returns a list of SymPy variables that describe the reactions that were chosen to be fixed when ensuring a
        steady-state solution exists. Note that fixed_reactions must be set to True in
        :func:`crnt4sbml.GeneralApproach.initialize_general_approach`.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> GA = network.get_general_approach()
        >>> signal = "C1"
        >>> response = "s1"
        >>> GA.initialize_general_approach(signal=signal, response=response, fix_reactions=True)
        >>> GA.get_fixed_reactions()
        """
        return self.__fixed_reactions

    def get_solutions_to_fixed_reactions(self):
        """
        Returns a list of SymPy expressions corresponding to the fixed reactions. The ordering of the elements
        corresponds to the list returned by :func:`crnt4sbml.GeneralApproach.get_fixed_reactions`. Note that
        fixed_reactions must be set to True in :func:`crnt4sbml.GeneralApproach.initialize_general_approach`.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> GA = network.get_general_approach()
        >>> signal = "C1"
        >>> response = "s1"
        >>> GA.initialize_general_approach(signal=signal, response=response, fix_reactions=True)
        >>> GA.get_solutions_to_fixed_reactions()
        """
        return self.__soln_to_fixed_reactions

    def get_conservation_laws(self):
        """
        Returns a string representation of the conservation laws. Here the values on the left hand side of each equation
        are the constants of the conservation laws.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> GA = network.get_general_approach()
        >>> print(GA.get_conservation_laws())
        """
        rhs = self.__cons_laws_sympy
        laws = ""
        for i in range(rhs.shape[0]):
            laws += 'C' + str(i + 1) + ' = ' + str(rhs[i]) + '\n'

        return laws

    def get_determinant_of_jacobian(self):
        """
        Returns a Sympy expression of the determinant of the Jacobian, where the Jacobian is with respect to the
        independent species.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> GA = network.get_general_approach()
        >>> signal = "C1"
        >>> response = "s1"
        >>> GA.initialize_general_approach(signal=signal, response=response)
        >>> GA.get_determinant_of_jacobian()
        """
        return self.__det_jac

    def get_jacobian(self):
        """
        Returns a Sympy expression of the Jacobian, where the Jacobian is with respect to the independent species.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> GA = network.get_general_approach()
        >>> signal = "C1"
        >>> response = "s1"
        >>> GA.initialize_general_approach(signal=signal, response=response)
        >>> GA.get_jacobian()
        """
        return self.__jac_mod_system

    def get_comm(self):
        """
        Returns a mpi4py communicator if it has been initialized and None otherwise.
        """
        return self.__comm

    def get_my_rank(self):
        """
        Returns the rank assigned by mpi4py if it is initialized, otherwise None will be returned.
        """
        return self.__my_rank

