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


class GeneralApproach:
    """
    Class for constructing a more general approach to bistability detection for systems with mass action kinetics.
    """
    def __init__(self, cgraph):
        """
        Initialization of GeneralApproach class.

        See also
        ---------
        crnt4sbml.CRNT.get_general_approach
        """
        self.__cgraph = cgraph

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

    def __obj_func_fixed(self, x, h_func, x_full, lambda_cons_laws_sympy, bounds, full_constraints):

        count0 = 0
        for i in range(len(x_full) - len(lambda_cons_laws_sympy)):
            if i not in self.__fixed_reaction_indices:
                x_full[i] = x[count0]
                count0 += 1

        count = 0
        for i in lambda_cons_laws_sympy:
            x_full[len(self.__x_bar) + count] = i(*tuple(x[self.__R-len(self.__fixed_reaction_indices):self.__R-len(self.__fixed_reaction_indices) + self.__N]))
            count += 1

        for i in range(len(self.__lambda_fixed_reactions)):
            x_full[self.__fixed_reaction_indices[i]] = self.__lambda_fixed_reactions[i](*tuple([x_full[i] for i in range(len(x_full)) if i not in self.__fixed_reaction_indices]))

        fun = h_func(*tuple(x_full))

        # bounds of the optimization problem
        sumval = numpy.float64(0.0)
        for i in range(len(x_full) - len(lambda_cons_laws_sympy)):
            sumval += numpy.maximum(numpy.float64(0.0), numpy.float64(bounds[i][0]) - x_full[i])
            sumval += numpy.maximum(numpy.float64(0.0), x_full[i] - numpy.float64(bounds[i][1]))

        fun += sumval

        # constraints of the optimization problem
        if full_constraints:
            for i in range(len(full_constraints)):
                fun += full_constraints[i][0](x_full[0:len(x_full) - len(lambda_cons_laws_sympy)], full_constraints[i][1])

        if numpy.isnan(fun):
            return numpy.Inf
        else:
            return fun

    def __obj_func(self, x, h_func, x_full, lambda_cons_laws_sympy, bounds, full_constraints):

        x_full[0:len(self.__x_bar)] = x

        count = 0
        for i in lambda_cons_laws_sympy:
            x_full[len(self.__x_bar) + count] = i(*tuple(x[self.__R:self.__R + self.__N]))
            count += 1

        fun = h_func(*tuple(x_full))

        # bounds of the optimization problem
        sumval = numpy.float64(0.0)
        for i in range(len(x)):
            sumval += numpy.maximum(numpy.float64(0.0), numpy.float64(bounds[i][0]) - x[i])
            sumval += numpy.maximum(numpy.float64(0.0), x[i] - numpy.float64(bounds[i][1]))

        fun += sumval

        # constraints of the optimization problem
        if full_constraints:
            for i in range(len(full_constraints)):
                fun += full_constraints[i][0](x_full[0:len(x_full) - len(lambda_cons_laws_sympy)], full_constraints[i][1])

        if numpy.isnan(fun):
            return numpy.Inf
        else:
            return fun

    def __gather_single_value(self, value, number_of_values):

        temp_full_value = self.__comm.gather(value, root=0)

        if self.__my_rank == 0:
            full_value = numpy.zeros(number_of_values, dtype=numpy.float64)

            # putting all the obtained minimum values into a single array
            count = 0
            for i in temp_full_value:
                for j in i:
                    full_value[count] = j
                    count += 1
        else:
            full_value = None

        return full_value

    def __gather_list_of_values(self, values):

        full_values = self.__comm.gather(values, root=0)

        if self.__my_rank == 0:
            list_of_values = []
            for i in range(len(full_values)):
                list_of_values += full_values[i]
        else:
            list_of_values = []

        return list_of_values

    def __gather_numpy_array_of_values(self, values):

        full_values = self.__comm.gather(values, root=0)

        if self.__my_rank == 0:
            list_of_values = []
            for i in range(len(full_values)):
                list_of_values += full_values[i]
            array_of_values = numpy.zeros((len(list_of_values), list_of_values[0].shape[0]), dtype=numpy.float64)
            for i in range(len(list_of_values)):
                array_of_values[i, :] = list_of_values[i]
        else:
            array_of_values = None

        return array_of_values

    def __distribute_points(self, samples):

        if self.__my_rank == 0:

            # number of tasks per core
            tasks = len(samples) // self.__num_cores  # // calculates the floor

            # remainder
            r = len(samples) - self.__num_cores * tasks

            # array that holds how many tasks each core has
            tasks_core = numpy.zeros(self.__num_cores, dtype=numpy.int64)
            tasks_core.fill(tasks)

            # distributing in the remainder
            ii = 0
            while r > 0:
                tasks_core[ii] += 1
                r -= 1
                ii += 1

            sample_portion = samples[0:tasks_core[0], :]

            if self.__num_cores > 1:
                for i in range(1, self.__num_cores):
                    start = sum(tasks_core[0:i])
                    end = start + tasks_core[i]
                    self.__comm.send(samples[start:end, :], dest=i, tag=i * 11)

        else:
            if self.__num_cores > 1:
                sample_portion = self.__comm.recv(source=0, tag=self.__my_rank * 11)

        return sample_portion

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
                The number of iterations to run the feasible point method.
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

        if parallel_flag:

            # initializing MPI proccess
            from mpi4py import MPI
            self.__comm = MPI.COMM_WORLD
            self.__my_rank = self.__comm.Get_rank()
            self.__num_cores = self.__comm.Get_size()

            self.__comm.Barrier()

            start_time = MPI.Wtime()

            if self.__my_rank == 0:
                print("Starting optimization ...")

            # creating initial decision vectors for feasible point method
            if self.__my_rank == 0:
                samples, temp_bounds, x_full, feasible_point_sets, full_constraints = self.__initialize_optimization(
                    seed, iterations, bounds, constraints)
            else:
                samples = None
                temp_bounds = None
                x_full = None
                full_constraints = None

            sample_portion = self.__distribute_points(samples)

            # broadcast important variables
            temp_bounds = self.__comm.bcast(temp_bounds, root=0)
            x_full = self.__comm.bcast(x_full, root=0)
            full_constraints = self.__comm.bcast(full_constraints, root=0)

            self.__comm.Barrier()

            det_point_sets, det_point_sets_fun, obtained_minimums, smallest_value = self.__main_optimization_routine(sample_portion, temp_bounds, bounds, x_full,
                                                                                                                     full_constraints, confidence_level_flag, seed,
                                                                                                                     dual_annealing_iters, print_flag, parallel_flag)

            self.__comm.Barrier()

            if confidence_level_flag:

                full_obtained_minimums = self.__gather_single_value(obtained_minimums, iterations)

                if self.__my_rank == 0:
                    self.__confidence_level(full_obtained_minimums, change_in_rel_error)

            else:
                smallest_values = self.__comm.gather(smallest_value, root=0)
                if self.__my_rank == 0:
                    min_value = min(smallest_values)
                    self.__important_info += "Smallest value achieved by objective function: " + str(min_value)

            list_det_point_sets = self.__gather_list_of_values(det_point_sets)
            list_det_point_sets_fun = self.__gather_list_of_values(det_point_sets_fun)

            self.__comm.Barrier()

            if self.__my_rank == 0:
                end_time = MPI.Wtime()
                elapsed = end_time - start_time
                self.__important_info += str(len(list_det_point_sets)) + " point(s) passed the optimization criteria. " + "\n"
                print(f"Elapsed time for optimization in seconds: {elapsed}")

            return list_det_point_sets, list_det_point_sets_fun

        else:
            print("Starting optimization ...")
            start_t = time.time()
            samples, temp_bounds, x_full, feasible_point_sets, full_constraints = self.__initialize_optimization(seed, iterations, bounds, constraints)

            det_point_sets, det_point_sets_fun, obtained_minimums, smallest_value = self.__main_optimization_routine(feasible_point_sets, temp_bounds, bounds, x_full, full_constraints,
                                                                                                                     confidence_level_flag, seed, dual_annealing_iters, print_flag, parallel_flag)

            end_t = time.time()
            elapsed = end_t - start_t
            print("Elapsed time for optimization in seconds: " + str(elapsed))

            if confidence_level_flag:
                self.__confidence_level(obtained_minimums, change_in_rel_error)
                self.__important_info += str(len(det_point_sets_fun)) + " point(s) passed the optimization criteria. " + "\n"
            else:
                self.__important_info += "Smallest value achieved by objective function: " + str(smallest_value) + "\n"
                self.__important_info += str(len(det_point_sets_fun)) + " point(s) passed the optimization criteria. " + "\n"

            return det_point_sets, det_point_sets_fun

    def __initialize_optimization(self, seed, iterations, bounds, constraints):

        numpy.random.seed(seed)

        if self.__fix_reactions:
            samples = numpy.random.rand(iterations, len(bounds) - len(self.__fixed_reaction_indices))

            temp_bounds = [bounds[i] for i in range(len(bounds)) if i not in self.__fixed_reaction_indices]

            ranges = numpy.asarray(temp_bounds, dtype=numpy.float64)
            samples = samples * (ranges[:, 1] - ranges[:, 0]) + ranges[:, 0]

        else:
            temp_bounds = None
            samples = numpy.random.rand(iterations, len(bounds))

            ranges = numpy.asarray(bounds, dtype=numpy.float64)
            samples = samples * (ranges[:, 1] - ranges[:, 0]) + ranges[:, 0]

        x_full = numpy.zeros(len(self.__lagrangian_vars), dtype=numpy.float64)

        feasible_point_sets = [samples[i] for i in range(iterations)]

        # setting up equality and inequality constraints provided by the user
        if constraints:
            full_constraints = []
            for i in constraints:
                if i["type"] == "ineq":
                    full_constraints.append(
                        [lambda x, func: numpy.maximum(numpy.float64(0.0), numpy.float64(-1.0 * func(x))), i["fun"]])          # TODO: make it be 0.5 times greater or something
                elif i["type"] == "eq":
                    full_constraints.append([lambda x, func: numpy.float64(numpy.abs(func(x))), i["fun"]])

                else:
                    print("The type of constraint provided is unknown. Please review the entered constraints.")
                    sys.exit()
        else:
            full_constraints = []

        return samples, temp_bounds, x_full, feasible_point_sets, full_constraints

    def __main_optimization_routine(self, feasible_point_sets, temp_bounds, bounds, x_full, full_constraints,
                                    confidence_level_flag, seed, dual_annealing_iters, print_flag, parallel_flag):

        det_point_sets = []
        det_point_sets_fun = []
        smallest_value = numpy.float(1e8)
        if confidence_level_flag:
            obtained_minimums = numpy.zeros(len(feasible_point_sets), dtype=numpy.float64)
        else:
            obtained_minimums = None

        if len(feasible_point_sets) != 0:

            for i in range(len(feasible_point_sets)):

                with numpy.errstate(divide='ignore', invalid='ignore'):

                    if self.__fix_reactions:
                        result = scipy.optimize.dual_annealing(self.__obj_func_fixed, bounds=temp_bounds,
                                                               args=(self.__lambda_obj_func_h, x_full,
                                                                     self.__cons_laws_sympy_lamb, bounds,
                                                                     full_constraints),
                                                               x0=feasible_point_sets[i], seed=seed,
                                                               local_search_options={'method': "Nelder-Mead"},
                                                               maxiter=dual_annealing_iters)
                    else:
                        result = scipy.optimize.dual_annealing(self.__obj_func, bounds=bounds,
                                                               args=(self.__lambda_obj_func_h, x_full,
                                                                     self.__cons_laws_sympy_lamb, bounds,
                                                                     full_constraints),
                                                               x0=feasible_point_sets[i], seed=seed,
                                                               local_search_options={'method': "Nelder-Mead"},
                                                               maxiter=dual_annealing_iters)

                    if print_flag:
                        print("Global function value: " + str(result.fun))
                        print("Decision vector used: ")
                        print(result.x)

                    if abs(result.fun) > numpy.float64(1e-100):

                        if self.__fix_reactions:
                            result1 = scipy.optimize.minimize(self.__obj_func_fixed, result.x, args=(self.__lambda_obj_func_h,
                                                                                               x_full,
                                                                                               self.__cons_laws_sympy_lamb,
                                                                                               bounds, full_constraints),
                                                              method='Nelder-Mead', tol=1e-16)
                        else:
                            result1 = scipy.optimize.minimize(self.__obj_func, result.x, args=(self.__lambda_obj_func_h,
                                                                                               x_full,
                                                                                               self.__cons_laws_sympy_lamb,
                                                                                               bounds,
                                                                                               full_constraints),
                                                              method='Nelder-Mead', tol=1e-16)

                        if print_flag:
                            print("Local function value: " + str(result1.fun))
                            print("Decision vector used: ")
                            print(result1.x)

                        if smallest_value > result1.fun:
                            smallest_value = result1.fun

                        if confidence_level_flag:
                            obtained_minimums[i] = result1.fun

                        if abs(result1.fun) <= numpy.finfo(float).eps:
                            # det_point_sets.append(result1.x)
                            out = self.__construct_full_reactions_and_species(result1.x, x_full, self.__cons_laws_sympy_lamb)
                            det_point_sets.append(out)
                            det_point_sets_fun.append(result1.fun)

                    else:
                        if smallest_value > result.fun:
                            smallest_value = result.fun
                        # det_point_sets.append(result.x)
                        out = self.__construct_full_reactions_and_species(result.x, x_full, self.__cons_laws_sympy_lamb)
                        det_point_sets.append(out)
                        det_point_sets_fun.append(result.fun)
                        if confidence_level_flag:
                            obtained_minimums[i] = result.fun

                    if print_flag:
                        print("")
        elif not parallel_flag:
            print("Determinant optimization has finished.")
            raise Exception("Optimization needs to be run with more iterations or different bounds.")

        return det_point_sets, det_point_sets_fun, obtained_minimums, smallest_value

    def __construct_full_reactions_and_species(self, x, x_full, lambda_cons_laws_sympy):

        if self.__fix_reactions:

            count0 = 0
            for i in range(len(x_full) - len(lambda_cons_laws_sympy)):
                if i not in self.__fixed_reaction_indices:
                    x_full[i] = x[count0]
                    count0 += 1

            count = 0
            for i in lambda_cons_laws_sympy:
                x_full[len(self.__x_bar) + count] = i(*tuple(x[self.__R - len(
                    self.__fixed_reaction_indices):self.__R - len(self.__fixed_reaction_indices) + self.__N]))
                count += 1

            for i in range(len(self.__lambda_fixed_reactions)):
                x_full[self.__fixed_reaction_indices[i]] = self.__lambda_fixed_reactions[i](
                    *tuple([x_full[i] for i in range(len(x_full)) if i not in self.__fixed_reaction_indices]))

            return numpy.array(list(x_full[0:len(self.__x_bar)]))

        else:
            return x

    def __confidence_level(self, obtained_minimums, change_in_rel_error):

        a = 1
        b = 5

        unique_elements, counts_elements = numpy.unique(obtained_minimums, return_counts=True)
        min_val_index = numpy.nanargmin(unique_elements)

        f_til = unique_elements[min_val_index]

        numpy_float64_smallest_positive_value = numpy.nextafter(numpy.float64(0), numpy.float64(1))

        if f_til > numpy_float64_smallest_positive_value:

            r = numpy.count_nonzero(
                numpy.abs(f_til - obtained_minimums) / f_til < numpy.float64(change_in_rel_error))

            n_til = obtained_minimums.shape[0]
            a_bar = a + b - 1
            b_bar = b - r - 1

            prob = 1 - (math.factorial(n_til + a_bar) * math.factorial(2 * n_til + b_bar)) / (
                    math.factorial(2 * n_til + a_bar) * math.factorial(n_til + b_bar))

        else:

            prob = 1.0

        self.__important_info += f"It was found that {unique_elements[min_val_index]} is the minimum objective function value with a confidence level of {prob} .\n"

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
                print("")
                print("A parallel version of numerical continuation is not available.")
                print("Please rerun your script without mpiexec.")
                print(
                    "For your convenience, the provided parameters have been saved in the current directory under the name params.npy.")
                numpy.save('./params.npy', parameters)

            sys.exit()

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

        self.__print_flag = print_lbls_flag

        multistable_param_ind, important_info, plot_specifications = BistabilityFinder.run_greedy_continuity_analysis \
            (species_num, parameters, self.__initialize_ant_string, self.__finalize_ant_string, species_y, dir_path,
             print_lbls_flag, auto_parameters, plot_labels)

        self.__important_info += important_info

        return multistable_param_ind, plot_specifications

    def __distribute_list_of_points(self, samples):

        if self.__my_rank == 0:

            # number of tasks per core
            tasks = len(samples) // self.__num_cores  # // calculates the floor

            # remainder
            r = len(samples) - self.__num_cores * tasks

            # array that holds how many tasks each core has
            tasks_core = numpy.zeros(self.__num_cores, dtype=numpy.int64)
            tasks_core.fill(tasks)

            # distributing in the remainder
            ii = 0
            while r > 0:
                tasks_core[ii] += 1
                r -= 1
                ii += 1

            sample_portion = samples[0:tasks_core[0]]

            if self.__num_cores > 1:
                for i in range(1, self.__num_cores):
                    start = sum(tasks_core[0:i])
                    end = start + tasks_core[i]
                    self.__comm.send(samples[start:end], dest=i, tag=i * 11)

        else:
            if self.__num_cores > 1:
                sample_portion = self.__comm.recv(source=0, tag=self.__my_rank * 11)

        return sample_portion

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
                print("")
                print("A parallel version of numerical continuation is not available.")
                print("Please rerun your script without mpiexec.")
                print(
                    "For your convenience, the provided parameters have been saved in the current directory under the name params.npy.")
                numpy.save('./params.npy', parameters)

            sys.exit()

        # setting default values for AUTO
        if 'NMX' not in auto_parameters.keys():
            auto_parameters['NMX'] = 10000

        if 'ITMX' not in auto_parameters.keys():
            auto_parameters['ITMX'] = 100

        # making the directory if it doesn't exist
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)

        self.__print_flag = print_lbls_flag

        species_num = [str(i) for i in self.__indp_species].index(species) + 1

        species_y = str(self.__indp_species[species_num - 1])

        multistable_param_ind, important_info, plot_specifications = BistabilityFinder.run_continuity_analysis \
            (species_num, parameters, self.__initialize_ant_string, self.__finalize_ant_string, species_y, dir_path,
             print_lbls_flag, auto_parameters, plot_labels)

        self.__important_info += important_info

        return multistable_param_ind, plot_specifications

    def __find_viable_indices(self, result_x, itg, spec_index, change_in_relative_error):

        conservation_vals = [self.__cons_laws_sympy_lamb[i](*tuple(result_x[self.__R:self.__R + self.__N]))
                             for i in range(len(self.__cons_laws_sympy_lamb))]

        cons_laws_spec_info = []
        for i in self.__cons_laws_sympy:
            spec_in_law = [ii for ii in self.__sympy_species if ii in i.atoms()]
            spec_indices = [self.__sympy_species.index(ii) for ii in spec_in_law]
            coeff_of_spec = [i.coeff(ii, 1) for ii in spec_in_law]
            cons_laws_spec_info.append([coeff_of_spec, spec_indices])

        temp_comb = list(itertools.product(*[i[1] for i in cons_laws_spec_info]))

        all_unique_comb = [temp_comb[i] for i in range(len(temp_comb)) if len(set(temp_comb[i])) == len(temp_comb[i])]

        viable_indices = []
        viable_out_values = []
        for i in all_unique_comb:
            initial_species_values = [0.0 for i in range(self.__N)]

            for j in range(len(conservation_vals)):
                initial_species_values[i[j]] = conservation_vals[j]

            out = self.__steady_state_finder(initial_species_values, result_x, spec_index, itg, change_in_relative_error)
            steady_cons = [self.__cons_laws_sympy_lamb[i](*tuple(out)) for i in range(len(self.__cons_laws_sympy_lamb))]

            if not numpy.array_equal(numpy.array(initial_species_values), out):

                # accepting those indices that are smaller than a predescribed relative error
                if all([abs(conservation_vals[ii] - steady_cons[ii])/abs(steady_cons[ii]) < change_in_relative_error for ii in range(len(conservation_vals))]):
                    viable_out_values.append(out)
                    viable_indices.append(i)
            elif len(out) == 2:

                # accepting those indices that are smaller than a predescribed relative error
                if all([abs(conservation_vals[ii] - steady_cons[ii]) / abs(steady_cons[ii]) < change_in_relative_error
                        for ii in range(len(conservation_vals))]):
                    viable_out_values.append(out)
                    viable_indices.append(i)

        return viable_indices, viable_out_values, conservation_vals

    def run_direct_simulation(self, params_for_global_min=None, dir_path="./", change_in_relative_error=1e-6, parallel_flag=False, print_flag=False, left_multiplier=0.5, right_multiplier=0.5):

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
        import scipy.integrate as itg

        lambda_inputs = self.__sympy_reactions + self.__sympy_species
        self.__ode_lambda_functions = [sympy.utilities.lambdify(lambda_inputs, self.__full_system[i]) for i in
                                       range(len(self.__full_system))]

        self.__jac_lambda_function = sympy.utilities.lambdify(lambda_inputs, self.__full_system.jacobian(self.__sympy_species))

        self.__dir_sim_print_flag = print_flag

        spec_index = self.__sympy_species.index(sympy.Symbol(self.__response, positive=True))

        if parallel_flag is False and self.__comm is None:

            #making the directory if it doesn't exist
            if not os.path.isdir(dir_path):
                os.mkdir(dir_path)

            print("Starting direct simulation ...")
            start_t = time.time()

        elif parallel_flag is True and self.__comm is None:

                from mpi4py import MPI
                self.__comm = MPI.COMM_WORLD
                self.__my_rank = self.__comm.Get_rank()
                self.__num_cores = self.__comm.Get_size()
                self.__comm.Barrier()

                if not os.path.isdir(dir_path) and self.__my_rank == 0:
                    os.mkdir(dir_path)

                self.__comm.Barrier()

                start_time = MPI.Wtime()

                if self.__my_rank == 0:
                    print("Starting direct simulation ...")
        elif self.__comm is not None:

            from mpi4py import MPI

            if not os.path.isdir(dir_path) and self.__my_rank == 0:
                os.mkdir(dir_path)

            self.__comm.Barrier()

            params_for_global_min = self.__comm.bcast(params_for_global_min, root=0)

            self.__comm.Barrier()

            start_time = MPI.Wtime()

            if self.__my_rank == 0:
                print("Starting direct simulation ...")

        else:
            print("Starting direct simulation ...")
            start_t = time.time()

        if len(params_for_global_min) == 0:
            print("The parameter sets provided has a length of zero, direct simulation cannot be ran.")
            sys.exit()

        viable_indices, viable_out_values, conservation_vals = self.__find_viable_indices(params_for_global_min[0], itg, spec_index, change_in_relative_error)

        spec_index, fwd_scan_vals, rvrs_scan_vals, fwd_scan_index, rvrs_scan_index = self.__initialize_direct_simulation(viable_indices, viable_out_values,
                                                                                                                         params_for_global_min[0], conservation_vals, itg,
                                                                                                                         change_in_relative_error, spec_index, left_multiplier,
                                                                                                                                 right_multiplier)
        if self.__dir_sim_print_flag:
            if self.__comm is None:
                self.__print_initial_conditions(fwd_scan_vals, rvrs_scan_vals)
            else:
                if self.__my_rank == 0:
                    self.__print_initial_conditions(fwd_scan_vals, rvrs_scan_vals)

        plot_flag = True

        for i in range(len(params_for_global_min)):

            if self.__dir_sim_print_flag:
                if self.__comm is None:
                    print(f"Conducting stability analysis of element {i} of the list provided ... ")
                else:
                    if self.__my_rank == 0:
                        print(f"Conducting stability analysis of element {i} of the list provided ... ")

            conservation_vals = [self.__cons_laws_sympy_lamb[ii](*tuple(params_for_global_min[i][self.__R:self.__R + self.__N]))
                                 for ii in range(len(self.__cons_laws_sympy_lamb))]

            con_law_value = conservation_vals[self.__signal_index]
            change_left = con_law_value * left_multiplier
            change_right = con_law_value * right_multiplier
            pcp_scan = numpy.linspace(con_law_value - change_left, con_law_value + change_right, 100)

            forward_scan, reverse_scan = self.__conduct_fwd_rvrs_scan(params_for_global_min[i], fwd_scan_vals,
                                                                      rvrs_scan_vals, pcp_scan, fwd_scan_index,
                                                                      rvrs_scan_index, spec_index, itg, change_in_relative_error)

            min_val, max_val = self.__get_min_max_vals(pcp_scan, forward_scan, reverse_scan)

            count = 0
            scan_vals = pcp_scan
            while count < 5:
                if len([ii for ii in scan_vals if ii >= min_val and ii <= max_val]) < 10:

                    second_scan = numpy.linspace(min_val, max_val, 60)

                    forward_scan, reverse_scan = self.__conduct_fwd_rvrs_scan(params_for_global_min[i], fwd_scan_vals,
                                                                              rvrs_scan_vals, second_scan, fwd_scan_index,
                                                                              rvrs_scan_index, spec_index, itg, change_in_relative_error)

                    min_val, max_val = self.__get_min_max_vals(second_scan, forward_scan, reverse_scan)
                    scan_vals = second_scan
                    count += 1
                else:
                    break

            if count == 0:
                second_scan = pcp_scan

            if plot_flag:

                if self.__comm is None:
                    self.__plot_direct_simulation(second_scan, forward_scan, reverse_scan, dir_path, i)
                else:
                    if self.__my_rank == 0:
                        self.__plot_direct_simulation(second_scan, forward_scan, reverse_scan, dir_path, i)


        if self.__comm is None:
            end_t = time.time()
            elapsed = end_t - start_t
            print("Elapsed time for direct simulation in seconds: " + str(elapsed))
        else:
            self.__comm.Barrier()
            if self.__my_rank == 0:
                end_time = MPI.Wtime()
                elapsed = end_time - start_time
                print(f"Elapsed time for direct simulation in seconds: {elapsed}")

    def __print_initial_conditions(self, fwd_scan_vals, rvrs_scan_vals):

        fwd_spec_inds = [i[0] for i in fwd_scan_vals]
        init_vals = []
        for i in range(self.__N):
            if i in fwd_spec_inds:
                init_vals.append(str(self.__sympy_species[i]) + " = " + "C" + str(fwd_spec_inds.index(i) + 1))
            else:
                init_vals.append(str(self.__sympy_species[i]) + " = 0.0")

        print(" ")
        print("For the forward scan the following initial condition will be used:")
        for i in init_vals:
            print(i)

        rvrs_spec_inds = [i[0] for i in rvrs_scan_vals]
        init_vals = []
        for i in range(self.__N):
            if i in rvrs_spec_inds:
                init_vals.append(str(self.__sympy_species[i]) + " = " + "C" + str(rvrs_spec_inds.index(i) + 1))
            else:
                init_vals.append(str(self.__sympy_species[i]) + " = 0.0")

        print(" ")
        print("For the reverse scan the following initial condition will be used:")
        for i in init_vals:
            print(i)
        print(" ")

    def __initialize_direct_simulation(self, viable_indices, viable_out_values, result_x, conservation_vals, itg,
                                       change_in_relative_error, spec_index, left_multiplier, right_multiplier):

        combos = list(itertools.combinations([i for i in range(len(viable_out_values))], 2))
        diff = [numpy.abs(viable_out_values[i[0]][spec_index] - viable_out_values[i[1]][spec_index]) for i in combos]
        maxpos = diff.index(max(diff))
        chosen_initial_combo = combos[maxpos]

        stop_flag = True
        while stop_flag:

            # selecting largest difference as the right pair
            maxpos = diff.index(max(diff))

            chosen_combo = combos[maxpos]

            fwd_scan_vals, rvrs_scan_vals, fwd_scan_index, rvrs_scan_index, fwd_ind, rvrs_ind = self.__get_important_scan_vals(chosen_combo,
                                                                                                                               spec_index,
                                                                                                                               viable_out_values,
                                                                                                                               viable_indices)

            # determining if the forward or reverse scan is constant, if so, remove it as a viable combination
            con_law_value = conservation_vals[self.__signal_index]
            # change = con_law_value*0.25
            # pcp_scan = numpy.linspace(con_law_value - change, con_law_value + change, 10)
            change_left = con_law_value * left_multiplier
            change_right = con_law_value * right_multiplier
            pcp_scan = numpy.linspace(con_law_value - change_left, con_law_value + change_right, 10)

            forward_scan, reverse_scan = self.__conduct_fwd_rvrs_scan(result_x, fwd_scan_vals,
                                                                      rvrs_scan_vals, pcp_scan, fwd_scan_index,
                                                                      rvrs_scan_index, spec_index, itg,
                                                                      change_in_relative_error)

            combos, diff, stop_flag, combos_flag = self.__get_new_combo(forward_scan, reverse_scan, fwd_ind, rvrs_ind, combos, diff)

            # if all combinations are thought to produce constant in time return the initial combo and continue
            if combos_flag:
                fwd_scan_vals, rvrs_scan_vals, fwd_scan_index, rvrs_scan_index, fwd_ind, rvrs_ind = self.__get_important_scan_vals(
                    chosen_initial_combo,
                    spec_index,
                    viable_out_values,
                    viable_indices)
                break

        return spec_index, fwd_scan_vals, rvrs_scan_vals, fwd_scan_index, rvrs_scan_index

    def __get_important_scan_vals(self, chosen_combo, spec_index, viable_out_values, viable_indices):

        # choosing the largest value at the species index as the "high concentration" option
        if viable_out_values[chosen_combo[0]][spec_index] < viable_out_values[chosen_combo[1]][spec_index]:

            fwd_scan_vals = [[viable_indices[chosen_combo[1]][j], j] for j in
                             range(len(viable_indices[chosen_combo[1]]))]
            fwd_ind = chosen_combo[1]
            rvrs_scan_vals = [[viable_indices[chosen_combo[0]][j], j] for j in
                              range(len(viable_indices[chosen_combo[0]]))]
            rvrs_ind = chosen_combo[0]

        else:
            fwd_scan_vals = [[viable_indices[chosen_combo[0]][j], j] for j in
                             range(len(viable_indices[chosen_combo[0]]))]
            fwd_ind = chosen_combo[0]
            rvrs_scan_vals = [[viable_indices[chosen_combo[1]][j], j] for j in
                              range(len(viable_indices[chosen_combo[1]]))]
            rvrs_ind = chosen_combo[1]

        # index to change in forward scan
        fwd_scan_index = [i[0] for i in fwd_scan_vals if i[1] == self.__signal_index][0]
        # index to change in reverse scan
        rvrs_scan_index = [i[0] for i in rvrs_scan_vals if i[1] == self.__signal_index][0]

        return fwd_scan_vals, rvrs_scan_vals, fwd_scan_index, rvrs_scan_index, fwd_ind, rvrs_ind

    def __get_new_combo(self, forward_scan, reverse_scan, fwd_ind, rvrs_ind, combos, diff):

        constant_index = []
        reverse_scan = [abs(i) for i in reverse_scan]
        forward_scan = [abs(i) for i in forward_scan]

        fwd_rel_change = abs(max(forward_scan) - min(forward_scan)) / max(forward_scan)

        if fwd_rel_change >= 0.98:
            constant_index.append(fwd_ind)

        rvrs_rel_change = abs(max(reverse_scan) - min(reverse_scan)) / max(reverse_scan)

        if rvrs_rel_change >= 0.98:
            constant_index.append(rvrs_ind)

        #if one or both of them are deemed to be constant, then we throw out one or both from the combos.
        #Use fwd_ind and rvrs_ind to throw out combos that were constant, then redo process with new
        #combo being created. Continue until process with no constant combo is found. If all combos are eleminated
        #return the initial combo.
        if constant_index:
            ind_to_remove = [i for i in range(len(combos)) if combos[i][0] in constant_index or combos[i][1] in constant_index]
            combos = [combos[i] for i in range(len(combos)) if i not in ind_to_remove]
            diff = [diff[i] for i in range(len(diff)) if i not in ind_to_remove]

            if not combos:
                return combos, diff, True, True

            return combos, diff, True, False
        else:
            return combos, diff, False, False

    def __get_min_max_vals(self, pcp_scan, forward_scan, reverse_scan):

        fwd_diff = [abs(forward_scan[i] - forward_scan[i + 1]) for i in range(len(forward_scan) - 1)]

        fwd_maxpos = fwd_diff.index(max(fwd_diff))

        if fwd_maxpos == 0:
            fwd_maxpos = 1
        elif fwd_maxpos+2 == len(pcp_scan):
            fwd_maxpos = len(pcp_scan)-2

        fwd_pcp_scan = pcp_scan[fwd_maxpos - 1:fwd_maxpos + 2]

        rvrs_diff = [abs(reverse_scan[i] - reverse_scan[i + 1]) for i in range(len(reverse_scan) - 1)]

        rvrs_maxpos = rvrs_diff.index(max(rvrs_diff))

        if rvrs_maxpos == 0:
            rvrs_maxpos = 1
        elif rvrs_maxpos+2 == len(pcp_scan):
            rvrs_maxpos = len(pcp_scan)-2

        rvrs_pcp_scan = pcp_scan[rvrs_maxpos - 1:rvrs_maxpos + 2]

        min_val = min(list(fwd_pcp_scan) + list(rvrs_pcp_scan))

        max_val = max(list(fwd_pcp_scan) + list(rvrs_pcp_scan))

        return min_val, max_val

    def __conduct_fwd_rvrs_scan(self, result_x, fwd_scan_vals, rvrs_scan_vals, pcp_scan, fwd_scan_index,
                                rvrs_scan_index, spec_index, itg, change_in_relative_error):

        if self.__comm is not None:
            pcp_scan = self.__distribute_list_of_points(pcp_scan)

        conservation_vals = [self.__cons_laws_sympy_lamb[i](*tuple(result_x[self.__R:self.__R + self.__N]))
                             for i in range(len(self.__cons_laws_sympy_lamb))]

        initial_species_values = [0.0 for i in range(self.__N)]

        for i in fwd_scan_vals:
            initial_species_values[i[0]] = conservation_vals[i[1]]

        forward_scan = []
        for i in pcp_scan:
            initial_species_values[fwd_scan_index] = i
            steady_state = self.__steady_state_finder(initial_species_values, result_x, spec_index, itg, change_in_relative_error)
            forward_scan.append(steady_state[spec_index])

        initial_species_values = [0.0 for i in range(self.__N)]
        for i in rvrs_scan_vals:
            initial_species_values[i[0]] = conservation_vals[i[1]]

        reverse_scan = []
        for i in pcp_scan:
            initial_species_values[rvrs_scan_index] = i
            steady_state = self.__steady_state_finder(initial_species_values, result_x, spec_index, itg, change_in_relative_error)
            reverse_scan.append(steady_state[spec_index])

        if self.__comm is not None:
            list_forward_scan = self.__gather_list_of_values(forward_scan)
            list_reverse_scan = self.__gather_list_of_values(reverse_scan)

            list_forward_scan = self.__comm.bcast(list_forward_scan, root=0)
            list_reverse_scan = self.__comm.bcast(list_reverse_scan, root=0)

            self.__comm.Barrier()

            return list_forward_scan, list_reverse_scan

        else:
            return forward_scan, reverse_scan

    def __steady_state_finder(self, initial_species_values, result_x, spec_index, itg, change_in_relative_error):

        def ff(t, cs, ks, ode_lambda_functions, jacobian):
            return [i(*tuple(ks), *tuple(cs)) for i in ode_lambda_functions]

        def jac_f(t, cs, ks, ode_lambda_functions, jacobian):
            return jacobian(*tuple(ks), *tuple(cs))

        len_time_interval = 100.0

        with numpy.errstate(divide='ignore', invalid='ignore'):
            out = itg.solve_ivp(ff, [0.0, len_time_interval], initial_species_values, args=(result_x[0:self.__R], self.__ode_lambda_functions, self.__jac_lambda_function), jac=jac_f, method='BDF', rtol=1e-6, atol=1e-9, vectorized=True) #'RK45')  #'LSODA')
            y0 = out.y[:, -1]

        flag = True

        i = 1
        while flag:
            tspan = [0.0 + i*len_time_interval, len_time_interval + i*len_time_interval]
            try:
                with numpy.errstate(divide='ignore', invalid='ignore'):
                    out = itg.solve_ivp(ff, tspan, y0, args=(result_x[0:self.__R], self.__ode_lambda_functions, self.__jac_lambda_function), jac=jac_f, method='BDF', rtol=1e-6, atol=1e-9, vectorized=True) #'RK45')  #'LSODA')
                    y0 = out.y[:, -1]
                i += 1
            except Exception as e:
                flag = False

            # if there is a division by zero we exit the routine
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    flag = abs(out.y[spec_index, -1] - out.y[spec_index, 0]) / abs(out.y[spec_index, -1]) > change_in_relative_error and i < 1000
            except Exception as e:
                flag = False

        return out.y[:, -1]

    def __plot_direct_simulation(self, pcp_scan, forward_scan, reverse_scan, path, index):

        import pandas
        from plotnine import ggplot, aes, geom_line, ylim, scale_color_distiller, facet_wrap, theme_bw, geom_path, \
            geom_point, labs, annotate
        from matplotlib import rc
        rc('text', usetex=True)

        out = pandas.DataFrame(columns=['dir', 'signal'] + [self.__response])
        for i in range(len(forward_scan)):
            out_i = pandas.DataFrame([forward_scan[i]], columns=[out.columns[2]])
            out_i['signal'] = pcp_scan[i]
            out_i['dir'] = 'Forward scan'
            out = pandas.concat([out, out_i[out.columns]])
        for i in range(len(reverse_scan)):
            out_i = pandas.DataFrame([reverse_scan[i]], columns=[out.columns[2]])
            out_i['signal'] = pcp_scan[i]
            out_i['dir'] = 'Reverse scan'
            out = pandas.concat([out, out_i[out.columns]])

        g = (ggplot(out)
             + aes(x='signal', y=self.__response, color='dir')
             + labs(x=f"{self.__signal} total", y=f"[{self.__response}]", color="")
             + geom_path(size=2, alpha=0.5)
             + geom_point(color="black")
             + theme_bw()
             + geom_point(color="black"))
        g.save(filename=path + f"/sim_bif_diag_{index}.png", format="png", width=6, height=4, units='in', verbose=False)


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

        if self.__print_flag:
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

    @staticmethod
    def __get_physiological_range(for_what=None):
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

                print("")
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

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> GA = network.get_general_approach()
        >>> signal = "C1"
        >>> response = "s1"
        >>> GA.initialize_general_approach(signal=signal, response=response)
        >>> GA.get_comm()
        """
        return self.__comm

    def get_my_rank(self):
        """
        Returns the rank assigned by mpi4py if it is initialized, otherwise None will be returned.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/sbml_file.xml")
        >>> GA = network.get_general_approach()
        >>> signal = "C1"
        >>> response = "s1"
        >>> GA.initialize_general_approach(signal=signal, response=response)
        >>> GA.get_my_rank()
        """
        return self.__my_rank

