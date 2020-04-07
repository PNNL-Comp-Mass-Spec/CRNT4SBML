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
        Initialization of GeneralApproach class using a network contructed by crnt4sbml.CRNT().

        See also
        ---------
        fill in this
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

    def initialize_general_approach(self, signal=None, response=None):

        if signal in ['C' + str(i + 1) for i in range(len(self.__cons_laws_sympy))] \
                and response in self.__cgraph.get_species():

            self.__signal = signal
            self.__response = response

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

    def __create_variables_for_optimization(self):

        self.__jac_mod_system = self.__indp_system_subs.jacobian(sympy.Matrix(self.__indp_species))

        # self.__jac_mod_system = self.__modified_ode_system.jacobian(sympy.Matrix(self.__species_mod_system))

        self.__det_jac = self.__jac_mod_system.det(method='lu')

        self.__build_lagrangian()

        self.__build_objective_function()

    def __build_lagrangian(self):

        self.__x_bar = self.__sympy_reactions + self.__sympy_species

        second_term = sympy.S.Zero
        for i in self.__indp_system_subs:
            second_term += i**2
        # for i in self.__modified_ode_system:
        #     second_term += i**2

        self.__lagrangian = self.__det_jac**2 + second_term

        self.__lagrangian_vars = self.__x_bar + [sympy.Symbol('C' + str(i + 1), real=True) for i in range(len(self.__cons_laws_sympy))]
        # self.__lagrangian_vars = self.__x_bar + [sympy.Symbol(self.__signal, real=True)]

        self.__lambda_lagrangian = sympy.utilities.lambdify(self.__lagrangian_vars, self.__lagrangian)

        # enforcing a steady state
        self.__steady_state_func = second_term

        self.__steady_state_lambda = sympy.utilities.lambdify(self.__lagrangian_vars, self.__steady_state_func)

    def __build_objective_function(self):

        self.__obj_func_h = sympy.S.Zero

        self.__obj_func_h = self.__lagrangian

        self.__lambda_obj_func_h = self.__lambda_lagrangian

    def __feasible_point_obj_func(self, x, steady_func, x_full, lambda_cons_laws_sympy, bounds):

        x_full[0:len(self.__x_bar)] = x

        # count = 0
        # for i in lambda_cons_laws_sympy:
        #     x_full[len(self.__x_bar) + count] = i(*tuple(x[self.__R:self.__R + self.__N]))
        #     count += 1

        x_full[len(self.__x_bar)] = lambda_cons_laws_sympy[self.__signal_index](*tuple(x[self.__R:self.__R + self.__N]))

        fun = steady_func(*tuple(x_full))

        sumval = numpy.float64(0.0)
        for i in range(len(x)):
            sumval += numpy.maximum(numpy.float64(0.0), numpy.float64(bounds[i][0]) - x[i])
            sumval += numpy.maximum(numpy.float64(0.0), x[i] - numpy.float64(bounds[i][1]))

        fun += sumval

        if numpy.isnan(fun):
            return numpy.Inf
        else:
            return fun


    def __obj_func(self, x, h_func, x_full, lambda_cons_laws_sympy, bounds):

        x_full[0:len(self.__x_bar)] = x

        # x_full[len(self.__x_bar)] = lambda_cons_laws_sympy[self.__signal_index](*tuple(x[self.__R:self.__R + self.__N]))

        count = 0
        for i in lambda_cons_laws_sympy:
            x_full[len(self.__x_bar) + count] = i(*tuple(x[self.__R:self.__R + self.__N]))
            count += 1

        fun = h_func(*tuple(x_full))

        sumval = numpy.float64(0.0)
        for i in range(len(x)):
            sumval += numpy.maximum(numpy.float64(0.0), numpy.float64(bounds[i][0]) - x[i])
            sumval += numpy.maximum(numpy.float64(0.0), x[i] - numpy.float64(bounds[i][1]))

        fun += sumval

        if numpy.isnan(fun):
            return numpy.Inf
        else:
            return fun

    def run_mpi_optimization(self, bounds=None, iterations=10, seed=0, print_flag=False, dual_annealing_iters=1000,
                             confidence_level_flag=False, change_in_rel_error=1e-2):

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
            dual_annealing_iters: integer
                The number of iterations that should be ran for dual annealing routine in optimization.
            confidence_level_flag: bool
                If True a confidence level for the objective function will be given.
            change_in_rel_error: float
                The maximum relative error that should be allowed to consider :math:`f_k` in the neighborhood
                of :math:`\widetilde{f}`.
        Returns
        --------
        params_for_global_min: list of numpy arrays
            A list of numpy arrays that correspond to the decision vectors of the problem.
        obj_fun_val_for_params: list of floats
            A list of objective function values produced by the corresponding decision vectors in params_for_global_min.
        my_rank: integer
            An integer corresponding to the rank assigned to the core within the routine.

        Examples
        ---------
        Fill in
        """

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
            numpy.random.seed(seed)
            samples = numpy.random.rand(iterations, len(self.__decision_vector))
            ranges = numpy.asarray(decision_vector_bounds, dtype=numpy.float64)
            samples = samples * (ranges[:, 1] - ranges[:, 0]) + ranges[:, 0]

            # import math
            # ranges = numpy.asarray([(math.log10(i[0]), math.log10(i[1])) for i in decision_vector_bounds],
            #                        dtype=numpy.float64)
            # samples = samples * (ranges[:, 1] - ranges[:, 0]) + ranges[:, 0]
            # samples = numpy.power(10, samples)
        else:
            samples = None

        sample_portion = self.__distribute_points(samples)


        self.__comm.Barrier()

        if self.__my_rank == 0:
            print("Feasible point method has finished.")

        if self.__my_rank != 0:
            self.__important_info += f"The number of feasible points used in determinant optimization by core {self.__my_rank}: " \
                                    + str(len(feasible_point_sets))
        else:
            self.__important_info += f"The number of feasible points used in determinant optimization by core {self.__my_rank}: " \
                                     + str(len(feasible_point_sets)) + "\n"

        det_point_sets, det_point_sets_fun, obtained_minimums, full_set_of_values, inputs, input_values, smallest_value = \
            self.__main_optimization_routine(decision_vector_bounds, bounds, feasible_point_sets, fixed_reaction_ind_all
                                             , indp_spec_ind_dec, confidence_level_flag, seed, dual_annealing_iters,
                                             print_flag)

        self.__comm.Barrier()

        if self.__my_rank == 0:
            end_time = MPI.Wtime()
            elapsed = end_time - start_time
            print(f"Elapsed time for optimization in seconds: {elapsed}")

        params = self.__final_check(det_point_sets, bounds, full_set_of_values, fixed_reaction_ind_all,
                                    self.__cons_laws_lamb, indp_spec_ind_dec, inputs, input_values)

        if confidence_level_flag:

            full_obtained_minimums = self.__gather_single_value(obtained_minimums, iterations)

            if self.__my_rank == 0:
                self.__confidence_level(full_obtained_minimums, change_in_rel_error)

        else:
            smallest_values = self.__comm.gather(smallest_value, root=0)
            if self.__my_rank == 0:
                min_value = min(smallest_values)
                self.__important_info += "Smallest value achieved by objective function: " + str(min_value)

        list_params = self.__gather_list_of_values(params)
        list_det_point_sets_fun = self.__gather_list_of_values(det_point_sets_fun)

        self.__comm.Barrier()

        return list_params, list_det_point_sets_fun, self.__my_rank

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
                         confidence_level_flag=False, change_in_rel_error=1e-2):

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
            dual_annealing_iters: integer
                The number of iterations that should be ran for dual annealing routine in optimization.
            confidence_level_flag: bool
                If True a confidence level for the objective function will be given.
            change_in_rel_error: float
                The maximum relative error that should be allowed to consider :math:`f_k` in the neighborhood
                of :math:`\widetilde{f}`.
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

        print("Running feasible point method for " + str(iterations) + " iterations ...")
        feasible_point_sets, x_full = self.__feasible_point_method(bounds, iterations, seed, print_flag, confidence_level_flag, dual_annealing_iters)
        print("Feasible point method has finished.")

        self.__important_info += "\nThe number of feasible points used in determinant optimization: " + str(len(feasible_point_sets)) + "\n"

        det_point_sets = []
        det_point_sets_fun = []
        smallest_value = numpy.float(1e8)
        if confidence_level_flag:
            obtained_minimums = numpy.zeros(iterations, dtype=numpy.float64)

        if len(feasible_point_sets) != 0:

            print("Starting determinant optimization ...")

            for i in range(len(feasible_point_sets)):

                with numpy.errstate(divide='ignore', invalid='ignore'):

                    result = scipy.optimize.dual_annealing(self.__obj_func, bounds=bounds, args=(self.__lambda_obj_func_h, x_full,
                                                                                                 self.__cons_laws_sympy_lamb, bounds), x0=feasible_point_sets[i], seed=seed,
                                                           local_search_options={'method': "Nelder-Mead"}, maxiter=dual_annealing_iters)

                    # result = scipy.optimize.basinhopping(self.__obj_func, feasible_point_sets[i],
                    #                                      minimizer_kwargs={'method': 'Nelder-Mead',
                    #                                                        'args': (self.__lambda_obj_func_h, x_full,
                    #                                                                 self.__cons_laws_sympy_lamb),
                    #                                                        'tol': 1e-16},
                    #                                      niter=dual_annealing_iters, seed=seed)

                    # result = scipy.optimize.differential_evolution(self.__obj_func, full_bounds,
                    #                                                args=(self.__lambda_obj_func_h, x_full))

                    if print_flag:
                        print("Global function value: " + str(result.fun))
                        print("Decision vector used: ")
                        print(result.x)

                    if abs(result.fun) > numpy.float64(1e-100):

                        result1 = scipy.optimize.minimize(self.__obj_func, result.x, args=(self.__lambda_obj_func_h,
                                                                                           x_full,
                                                                                           self.__cons_laws_sympy_lamb,
                                                                                           bounds),
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
                            det_point_sets.append(result1.x)
                            det_point_sets_fun.append(result1.fun)

                    else:
                        if smallest_value > result.fun:
                            smallest_value = result.fun
                        det_point_sets.append(result.x)
                        det_point_sets_fun.append(result.fun)
                        if confidence_level_flag:
                            obtained_minimums[i] = result.fun

                    if print_flag:
                        print("")
        else:
            print("Determinant optimization has finished.")
            raise Exception("Optimization needs to be run with more iterations or different bounds.")

        print("Determinant optimization has finished.")

        end_t = time.time()
        elapsed = end_t - start_t
        print("Elapsed time for optimization in seconds: " + str(elapsed))

        if confidence_level_flag:
            self.__confidence_level(obtained_minimums, change_in_rel_error)

        else:
            self.__important_info += "Smallest value achieved by objective function: " + str(smallest_value) + "\n"

        return det_point_sets, det_point_sets_fun

    def __feasible_point_method(self, bounds, iterations, seed, print_flag, confidence_level_flag, dual_annealing_iters):

        numpy.random.seed(seed)

        samples = numpy.random.rand(iterations, len(bounds))

        ranges = numpy.asarray(bounds, dtype=numpy.float64)
        samples = samples * (ranges[:, 1] - ranges[:, 0]) + ranges[:, 0]

        lower_bounds = numpy.zeros(len(bounds), dtype=numpy.float64)
        upper_bounds = numpy.zeros(len(bounds), dtype=numpy.float64)

        for i in range(len(bounds)):
            lower_bounds[i] = numpy.float64(bounds[i][0])
            upper_bounds[i] = numpy.float64(bounds[i][1])

        x_full = numpy.zeros(len(self.__lagrangian_vars), dtype=numpy.float64)

        # x_full[len(bounds) + 1: len(bounds) + len(bounds) + 1] = lower_bounds
        #
        # x_full[len(bounds) + len(bounds) + 1:] = upper_bounds

        # x_full[len(bounds) + len(self.__cons_laws_sympy): len(bounds) + len(bounds) + len(
        #     self.__cons_laws_sympy)] = lower_bounds
        #
        # x_full[len(bounds) + len(bounds) + len(self.__cons_laws_sympy):] = upper_bounds

        feasible_point_sets = []
        if False:
            for i in range(iterations):

                with numpy.errstate(divide='ignore', invalid='ignore'):

                    result = scipy.optimize.dual_annealing(self.__feasible_point_obj_func,
                                                           bounds=bounds,
                                                           args=(self.__steady_state_lambda, x_full, self.__cons_laws_sympy_lamb,
                                                                 bounds),
                                                           x0=samples[i],
                                                           seed=seed,
                                                           local_search_options={'method': "Nelder-Mead"}, maxiter=10) #dual_annealing_iters)

                    # result = scipy.optimize.basinhopping(self.__feasible_point_obj_func, samples[i],
                    #                                      minimizer_kwargs={'method': 'Nelder-Mead',
                    #                                                        'args': (self.__steady_state_lambda, x_full,
                    #                                                                 self.__cons_laws_sympy_lamb),
                    #                                                        'tol': 1e-16},
                    #                                      niter=10, seed=seed)

                    if print_flag:
                        print("Objective function value: " + str(result.fun))
                        print("Decision vector used: ")
                        print(result.x)

                    if abs(result.fun) > numpy.float64(1e-100):
                        result1 = scipy.optimize.minimize(self.__feasible_point_obj_func, result.x,
                                                          args=(self.__steady_state_lambda, x_full, self.__cons_laws_sympy_lamb,
                                                                bounds),
                                                          method='Nelder-Mead', tol=1e-16)

                        if print_flag:
                            print("Objective function value: " + str(result1.fun))
                            print("Decision vector used: ")
                            print(result1.x)

                        if abs(result1.fun) <= numpy.finfo(float).eps or confidence_level_flag:
                            feasible_point_sets.append(result1.x)

                    else:
                        feasible_point_sets.append(result.x)

                    if print_flag:
                        print("")
        feasible_point_sets = [samples[i] for i in range(iterations)]
        return feasible_point_sets, x_full




    def __main_optimization_routine(self, decision_vector_bounds, bounds, feasible_point_sets, fixed_reaction_ind_all,
                                    indp_spec_ind_dec, confidence_level_flag, seed, dual_annealing_iters, print_flag):

        self.__comm.Barrier()

        if self.__my_rank == 0:
            print("Starting determinant optimization ...")

        bounds_2 = decision_vector_bounds
        if self.__fix_reactions:
            inputs = numpy.zeros(len(self.__vars_for_lam_fixed_reactions), dtype=numpy.float64)
        else:
            inputs = numpy.zeros(len(self.__lambda_variables), dtype=numpy.float64)

        full_set_of_values = numpy.zeros(len(bounds), dtype=numpy.float64)

        input_values = numpy.zeros(len(self.__lambda_variables), dtype=numpy.float64)

        det_point_sets = []
        det_point_sets_fun = []
        smallest_value = numpy.float(1e8)

        obtained_minimums = numpy.zeros(len(feasible_point_sets), dtype=numpy.float64)

        if len(feasible_point_sets) != 0:

            for i in range(len(feasible_point_sets)):

                with numpy.errstate(divide='ignore', invalid='ignore'):

                    result = scipy.optimize.dual_annealing(self.__obj_func, bounds=bounds_2,
                                                           args=(self.__lambda_det_jac, bounds,
                                                                 full_set_of_values, fixed_reaction_ind_all,
                                                                 self.__cons_laws_lamb, indp_spec_ind_dec,
                                                                 inputs, input_values), x0=feasible_point_sets[i],
                                                           seed=seed,
                                                           local_search_options={'method': "Nelder-Mead"},
                                                           maxiter=dual_annealing_iters)

                    if print_flag:
                        print("Global function value: " + str(result.fun))
                        print("Decision vector used: ")
                        print(result.x)

                    if abs(result.fun) > numpy.float64(1e-100):

                        result1 = scipy.optimize.minimize(self.__obj_func, result.x,
                                                          args=(self.__lambda_det_jac, bounds,
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

                        if confidence_level_flag:
                            obtained_minimums[i] = result1.fun

                        if abs(result1.fun) <= numpy.finfo(float).eps:
                            det_point_sets.append(result1.x)
                            det_point_sets_fun.append(result1.fun)

                    else:
                        if smallest_value > result.fun:
                            smallest_value = result.fun
                        det_point_sets.append(result.x)
                        det_point_sets_fun.append(result.fun)
                        if confidence_level_flag:
                            obtained_minimums[i] = result.fun

                    if print_flag:
                        print("")

            self.__comm.Barrier()
            if self.__my_rank == 0:
                print("Determinant optimization has finished.")


            return det_point_sets, det_point_sets_fun, obtained_minimums, full_set_of_values, inputs, input_values, smallest_value

        else:
            raise Exception("Optimization needs to be run with more iterations or different bounds.")


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
        Function for running the greedy numerical continuation and bistability analysis portions of the mass conservation
        approach. This routine uses the initial value of the principal continuation parameter to construct AUTO
        parameters and then tests varying fixed step sizes for the continuation problem. Note that this routine may
        produce jagged or missing sections in the plots provided. To produce better plots one should use the information
        provided by this routine to run :func:`crnt4sbml.GeneralApproach.run_continuity_analysis`.

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
        See .
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

        multistable_param_ind, important_info, plot_specifications = BistabilityFinder.run_greedy_continuity_analysis \
            (species_num, parameters, self.__initialize_ant_string, self.__finalize_ant_string, species_y, dir_path,
             print_lbls_flag, auto_parameters, plot_labels)

        self.__important_info += important_info

        return multistable_param_ind, plot_specifications

    # def run_mpi_greedy_continuity_analysis(self, species=None, parameters=None, dir_path="./num_cont_graphs",
    #                                        print_lbls_flag=False, auto_parameters=None, plot_labels=None):
    #
    #     """
    #     Function for running the greedy numerical continuation and bistability analysis portions of the mass conservation
    #     approach. This routine uses the initial value of the principal continuation parameter to construct AUTO
    #     parameters and then tests varying fixed step sizes for the continuation problem. Note that this routine may
    #     produce jagged or missing sections in the plots provided. To produce better plots one should use the information
    #     provided by this routine to run :func:`crnt4sbml.GeneralApproach.run_continuity_analysis`.
    #
    #     Parameters
    #     ------------
    #         species: string
    #             A string stating the species that is the y-axis of the bifurcation diagram.
    #         parameters: list of numpy arrays
    #             A list of numpy arrays corresponding to the decision vectors that produce a small objective function
    #             value.
    #         dir_path: string
    #             A string stating the path where the bifurcation diagrams should be saved.
    #         print_lbls_flag: bool
    #             If True the routine will print the special points found by AUTO 2000 and False will not print any
    #             special points.
    #         auto_parameters: dict
    #             Dictionary defining the parameters for the AUTO 2000 run. Please note that only the
    #             PrincipalContinuationParameter in this dictionary should be defined, no other AUTO parameters should
    #             be set. For more information on these parameters refer to :download:`AUTO parameters <../auto2000_input.pdf>`.
    #         plot_labels: list of strings
    #             A list of strings defining the labels for the x-axis, y-axis, and title. Where the first element
    #             is the label for x-axis, second is the y-axis label, and the last element is the title label. If
    #             you would like to use the default settings for some of the labels, simply provide None for that
    #             element.
    #     Returns
    #     ---------
    #         multistable_param_ind: list of integers
    #             A list of those indices in 'parameters' that produce multistable plots.
    #         sample_portion: list of 1D numpy arrays
    #             A list of 1D numpy arrays corresponding to those values in the input variable parameters that was
    #             distributed to the core.
    #         plot_specifications: list of lists
    #             A list whose elements correspond to the plot specifications of each element in multistable_param_ind.
    #             Each element is a list where the first element specifies the range used for the x-axis, the second
    #             element is the range for the y-axis, and the last element provides the x-y values and special point label
    #             for each special point in the plot.
    #
    #     Example
    #     ---------
    #     See .
    #     """
    #
    #     # setting default values for AUTO
    #     if 'NMX' not in auto_parameters.keys():
    #         auto_parameters['NMX'] = 10000
    #
    #     if 'ITMX' not in auto_parameters.keys():
    #         auto_parameters['ITMX'] = 100
    #
    #     # making the directory if it doesn't exist
    #     if not os.path.isdir(dir_path) and self.__my_rank == 0:
    #         os.mkdir(dir_path)
    #
    #     species_num = [str(i) for i in self.__indp_species].index(species) + 1
    #
    #     species_y = str(self.__indp_species[species_num-1])
    #
    #     if self.__comm is not None:
    #         sample_portion = self.__distribute_list_of_points(parameters)
    #         self.__comm.Barrier()
    #     else:
    #
    #         from mpi4py import MPI
    #         self.__comm = MPI.COMM_WORLD
    #         self.__my_rank = self.__comm.Get_rank()
    #         self.__num_cores = self.__comm.Get_size()
    #         self.__comm.Barrier()
    #
    #         if not os.path.isdir(dir_path) and self.__my_rank == 0:
    #             os.mkdir(dir_path)
    #         sample_portion = self.__distribute_list_of_points(parameters)
    #         self.__comm.Barrier()
    #
    #     multistable_param_ind, important_info, plot_specifications = BistabilityFinder.run_mpi_greedy_continuity_analysis \
    #         (species_num, sample_portion, self.__initialize_ant_string, self.__finalize_ant_string, species_y, dir_path,
    #          print_lbls_flag, auto_parameters, plot_labels, self.__my_rank, self.__comm)
    #
    #     self.__important_info += important_info
    #
    #     return multistable_param_ind, sample_portion, plot_specifications

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

    # def run_mpi_continuity_analysis(self, species=None, parameters=None, dir_path="./num_cont_graphs",
    #                                 print_lbls_flag=False, auto_parameters=None, plot_labels=None):
    #
    #     """
    #     Function for running the numerical continuation and bistability analysis portions of the mass conservation
    #     approach.
    #
    #     Parameters
    #     ------------
    #         species: string
    #             A string stating the species that is the y-axis of the bifurcation diagram.
    #         parameters: list of numpy arrays
    #             A list of numpy arrays corresponding to the decision vectors that produce a small objective function
    #             value.
    #         dir_path: string
    #             A string stating the path where the bifurcation diagrams should be saved.
    #         print_lbls_flag: bool
    #             If True the routine will print the special points found by AUTO 2000 and False will not print any
    #             special points.
    #         auto_parameters: dict
    #             Dictionary defining the parameters for the AUTO 2000 run. Please note that one should **not** set
    #             'SBML' or 'ScanDirection' in these parameters as these are automatically assigned. It is absolutely
    #             necessary to set PrincipalContinuationParameter in this dictionary. For more information on these
    #             parameters refer to :download:`AUTO parameters <../auto2000_input.pdf>`. 'NMX' will default to
    #             10000 and 'ITMX' to 100.
    #         plot_labels: list of strings
    #             A list of strings defining the labels for the x-axis, y-axis, and title. Where the first element
    #             is the label for x-axis, second is the y-axis label, and the last element is the title label. If
    #             you would like to use the default settings for some of the labels, simply provide None for that
    #             element.
    #     Returns
    #     ---------
    #         multistable_param_ind: list of integers
    #             A list of those indices in 'parameters' that produce multistable plots.
    #         sample_portion: list of 1D numpy arrays
    #             A list of 1D numpy arrays corresponding to those values in the input variable parameters that was
    #             distributed to the core.
    #         plot_specifications: list of lists
    #             A list whose elements correspond to the plot specifications of each element in multistable_param_ind.
    #             Each element is a list where the first element specifies the range used for the x-axis, the second
    #             element is the range for the y-axis, and the last element provides the x-y values and special point label
    #             for each special point in the plot.
    #
    #     Example
    #     ---------
    #     See .
    #     """
    #
    #     # setting default values for AUTO
    #     if 'NMX' not in auto_parameters.keys():
    #         auto_parameters['NMX'] = 10000
    #
    #     if 'ITMX' not in auto_parameters.keys():
    #         auto_parameters['ITMX'] = 100
    #
    #     # making the directory if it doesn't exist
    #     if not os.path.isdir(dir_path) and self.__my_rank == 0:
    #         os.mkdir(dir_path)
    #
    #     species_num = [str(i) for i in self.__indp_species].index(species) + 1
    #
    #     species_y = str(self.__indp_species[species_num - 1])
    #
    #     if self.__comm is not None:
    #         sample_portion = self.__distribute_list_of_points(parameters)
    #         self.__comm.Barrier()
    #     else:
    #
    #         from mpi4py import MPI
    #         self.__comm = MPI.COMM_WORLD
    #         self.__my_rank = self.__comm.Get_rank()
    #         self.__num_cores = self.__comm.Get_size()
    #         self.__comm.Barrier()
    #
    #         if not os.path.isdir(dir_path) and self.__my_rank == 0:
    #             os.mkdir(dir_path)
    #         sample_portion = self.__distribute_list_of_points(parameters)
    #         self.__comm.Barrier()
    #
    #     multistable_param_ind, important_info, plot_specifications = BistabilityFinder.run_mpi_continuity_analysis \
    #         (species_num, sample_portion, self.__initialize_ant_string, self.__finalize_ant_string, species_y, dir_path,
    #          print_lbls_flag, auto_parameters, plot_labels, self.__my_rank, self.__comm)
    #
    #     self.__important_info += important_info
    #
    #     return multistable_param_ind, sample_portion, plot_specifications

    def run_continuity_analysis(self, species=None, parameters=None, dir_path="./num_cont_graphs",
                                       print_lbls_flag=False, auto_parameters=None, plot_labels=None):

        """
        Function for running the numerical continuation and bistability analysis portions of the mass conservation
        approach.

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
        See .
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

        species_y = str(self.__indp_species[species_num - 1])

        multistable_param_ind, important_info, plot_specifications = BistabilityFinder.run_continuity_analysis \
            (species_num, parameters, self.__initialize_ant_string, self.__finalize_ant_string, species_y, dir_path,
             print_lbls_flag, auto_parameters, plot_labels)

        self.__important_info += important_info

        return multistable_param_ind, plot_specifications

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
        print(ant_str)

        return ant_str

    def get_optimization_bounds(self):
        """
        Returns a list of tuples that corresponds to the determined physiological bounds chosen for the problem. Each
        entry corresponds to the list provided by :func:`crnt4sbml.GeneralApproach.get_input_vector`.

         Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_optimization_bounds()
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

    def generate_report(self):
        """
        Prints out helpful details constructed by :func:`crnt4sbml.GeneralApproach.run_optimization` and
        :func:`crnt4sbml.GeneralApproach.run_continuity_analysis`.

        Example
        --------
        See .
        """
        if self.__comm == None:
            print(self.__important_info)
        else:

            all_important_info = self.__comm.gather(self.__important_info, root=0)
            self.__comm.Barrier()

            if self.__my_rank == 0:

                print("")

                for i in range(1, len(all_important_info)):
                    print(all_important_info[i])
                print(self.__important_info)

    def get_input_vector(self):
        """
        Returns a list of SymPy variables that specifies the ordering of the reactions and species for which bounds
        need to be provided.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_input_vector()
        """
        return self.__sympy_reactions + self.__sympy_species

    def get_decision_vector(self):
        """
        Returns a list of SymPy variables that specifies the ordering of the reactions and species of the decision
        vector to be used in optimization. Note: this method should not be used to create bounds for the optimization
        routine, rather :func:`crnt4sbml.GeneralApproach.get_input_vector` should be used.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_decision_vector()
        """
        return self.__decision_vector

    def get_independent_odes_subs(self):
        """
        Returns a Sympy Matrix representing the independent ODE system with conservation laws substituted in. Each row
        corresponds to the ODE for the species corresponding to the list provided by :func:`crnt4sbml.GeneralApproach.get_independent_species`.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_independent_odes_subs()
        """
        return self.__indp_system_subs

    def get_independent_odes(self):
        """
        Returns a Sympy Matrix representing the independent ODE system without conservation laws substituted in. Each row
        corresponds to the ODE for the species corresponding to the list provided by :func:`crnt4sbml.GeneralApproach.get_independent_species`.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_independent_odes()
        """
        return self.__indp_system

    def get_independent_species(self):
        """
        Returns a list of SymPy variables that reflects the independent species chosen for the general approach.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_independent_species()
        """

        return self.__indp_species

    def get_fixed_reactions(self):
        """
        Returns a list of SymPy variables that describe the reactions that were chosen to be fixed when ensuring a
        steady-state solution exists.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_fixed_reactions()
        """
        return self.__fixed_reactions

    def get_solutions_to_fixed_reactions(self):
        """
        Returns a list of SymPy expressions corresponding to the fixed reactions. The ordering of the elements
        corresponds to the list returned by :func:`crnt4sbml.GeneralApproach.get_fixed_reactions`.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_solutions_to_fixed_reactions()
        """
        return self.__soln_to_fixed_reactions2

    def get_conservation_laws(self):
        """
        Returns a string representation of the conservation laws. Here the values on the left hand side of each equation
        are the constants of the conservation laws.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> print(approach.get_conservation_laws())
        """
        rhs = self.__cons_laws_sympy # self.__BT_mat * sympy.Matrix([self.__sympy_species]).T
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
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_determinant_of_jacobian()
        """
        return self.__det_jac_indp_subs

    def get_variables_for_lambda_functions(self):
        """
        Returns a list of all SymPy variables that are used in the independent ODE system with conservation laws
        substituted in, which is given by :func:`crnt4sbml.GeneralApproach.get_independent_odes_subs`.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_variables_for_lambda_functions()
        """
        return self.__lambda_variables

    def get_indpendent_odes_lambda_function(self):
        """
        Returns a list of lamda functions where each element corresponds to the row of the independent ODE system with
        conservation laws substituted in, which is given by :func:`crnt4sbml.GeneralApproach.get_independent_odes_subs`.
        The variables for the lambda functions are given by :func:`crnt4sbml.GeneralApproach.get_variables_for_lambda_functions`.

        Example
        --------
        >>> import crnt4sbml
        >>> network = crnt4sbml.CRNT("path/to/SBML_File.xml")
        >>> approach = network.get_general_approach(signal, response)
        >>> approach.get_indpendent_odes_lambda_function()
        """

        return self.__lambda_indp_odes

    def get_full_set_of_values(self, params_for_global_min):

        print(self.__lambda_variables)

    def get_lagrange_vars(self):

        return self.__lagrangian_vars

    def jac_modified_sys(self):

        return self.__jac_mod_system

    def lambda_jac_modified_sys(self):

        return sympy.utilities.lambdify(self.__lagrangian_vars, self.__jac_mod_system)



