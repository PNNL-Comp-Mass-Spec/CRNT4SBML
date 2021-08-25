# Parent class to MassConservationApproach, SemiDiffusiveApproach, and GeneralApproach
import scipy
import numpy.linalg
import time
import math
import sys


class BistabilityFinder(object):

    def __parent_run_optimization(self):

        # get expression needed to evaluate the current subclass variable
        self.__mangled_name = "self._" + self.__class__.__name__
        self.__method = eval(self.__mangled_name + "__method")
        self.__parallel_flag = eval(self.__mangled_name + "__parallel_flag")
        self.__confidence_level_flag = eval(self.__mangled_name + "__confidence_level_flag")
        self.__print_flag = eval(self.__mangled_name + "__print_flag")
        self.__change_in_rel_error = eval(self.__mangled_name + "__change_in_rel_error")
        self.__iterations = eval(self.__mangled_name + "__iterations")

        if self.__parallel_flag:
            samples = self.__initialize_mpi_optimization()

            det_point_sets, det_point_sets_fun, obtained_minimums, smallest_value = self.__main_optimization_routine(samples)

            return self.__finalize_mpi_optimization(det_point_sets, det_point_sets_fun, obtained_minimums, smallest_value)

        else:
            self.__my_rank = None
            self.__comm = None
            samples = self.__initialize_optimization()

            det_point_sets, det_point_sets_fun, obtained_minimums, smallest_value = self.__main_optimization_routine(samples)

            return self.__finalize_optimization(det_point_sets, det_point_sets_fun, obtained_minimums, smallest_value)


    def __initialize_optimization(self):

        # creating initial decision vectors for feasible point method
        if self.__method == "GeneralApproach":
            samples = self.__initialize_optimization_ga()
        elif self.__method == "MassConservationApproach":
            samples = self.__initialize_optimization_mca()
        elif self.__method == "SemiDiffusiveApproach":
            samples = self.__initialize_optimization_sda()

        return samples

    def __finalize_optimization(self, det_point_sets, det_point_sets_fun, obtained_minimums, smallest_value):

        important_info = ''
        if self.__confidence_level_flag:
            important_info = self.__confidence_level(obtained_minimums, self.__change_in_rel_error, important_info)
            important_info += str(len(det_point_sets_fun)) + " point(s) passed the optimization criteria. " + "\n"
        else:
            important_info += "Smallest value achieved by objective function: " + str(smallest_value) + "\n"
            important_info += str(len(det_point_sets_fun)) + " point(s) passed the optimization criteria. " + "\n"

        return det_point_sets, det_point_sets_fun, important_info

    def __initialize_mpi_optimization(self):

        # initializing MPI proccess
        global MPI
        from mpi4py import MPI
        global mpi_mod
        from .mpi_routines import MPIRoutines as mpi_mod

        self.__comm = MPI.COMM_WORLD
        self.__my_rank = self.__comm.Get_rank()
        self.__num_cores = self.__comm.Get_size()

        self.__comm.Barrier()

        # creating initial decision vectors for feasible point method
        if self.__method == "GeneralApproach":
            samples = self.__initialize_optimization_ga()
        elif self.__method == "MassConservationApproach":
            samples = self.__initialize_optimization_mca()
        elif self.__method == "SemiDiffusiveApproach":
            samples = self.__initialize_optimization_sda()

        return samples

    def __finalize_mpi_optimization(self, det_point_sets, det_point_sets_fun, obtained_minimums, smallest_value):

        important_info = ''

        self.__comm.Barrier()

        if self.__confidence_level_flag:

            full_obtained_minimums = mpi_mod.gather_single_value(obtained_minimums, self.__iterations, self.__comm, self.__my_rank)

            if self.__my_rank == 0:
                important_info = self.__confidence_level(full_obtained_minimums, self.__change_in_rel_error, important_info)

        else:
            smallest_values = self.__comm.gather(smallest_value, root=0)
            if self.__my_rank == 0:
                min_value = min(smallest_values)
                important_info += "Smallest value achieved by objective function: " + str(min_value) + "\n"

        list_det_point_sets = mpi_mod.gather_list_of_values(det_point_sets, self.__comm, self.__my_rank)
        list_det_point_sets_fun = mpi_mod.gather_list_of_values(det_point_sets_fun, self.__comm, self.__my_rank)

        self.__comm.Barrier()

        if self.__my_rank == 0:
            important_info += str(len(list_det_point_sets)) + " point(s) passed the optimization criteria. " + "\n"

        return list_det_point_sets, list_det_point_sets_fun, important_info

    def __initialize_optimization_ga(self):

        if self.__my_rank == 0 or self.__my_rank is None:

            numpy.random.seed(self._GeneralApproach__seed)

            if self._GeneralApproach__fix_reactions:
                samples = numpy.random.rand(self.__iterations, len(self._GeneralApproach__bounds) -
                                            len(self._GeneralApproach__fixed_reaction_indices))

                self._GeneralApproach__temp_bounds = [self._GeneralApproach__bounds[i] for i in
                                                      range(len(self._GeneralApproach__bounds)) if i not in
                                                      self._GeneralApproach__fixed_reaction_indices]

                ranges = numpy.asarray(self._GeneralApproach__temp_bounds, dtype=numpy.float64)
                samples = samples * (ranges[:, 1] - ranges[:, 0]) + ranges[:, 0]

            else:
                self._GeneralApproach__temp_bounds = None
                samples = numpy.random.rand(self.__iterations, len(self._GeneralApproach__bounds))

                ranges = numpy.asarray(self._GeneralApproach__bounds, dtype=numpy.float64)
                samples = samples * (ranges[:, 1] - ranges[:, 0]) + ranges[:, 0]

            self._GeneralApproach__x_full = numpy.zeros(len(self._GeneralApproach__lagrangian_vars), dtype=numpy.float64)

            # setting up equality and inequality constraints provided by the user
            if self._GeneralApproach__constraints:
                self._GeneralApproach__full_constraints = []
                for i in self._GeneralApproach__constraints:
                    if i["type"] == "ineq":
                        self._GeneralApproach__full_constraints.append([lambda x, func: numpy.maximum(numpy.float64(0.0), numpy.float64(-1.0 * func(x))),
                                                                        i["fun"]])   # TODO: make it be 0.5 times greater or something
                    elif i["type"] == "eq":
                        self._GeneralApproach__full_constraints.append([lambda x, func: numpy.float64(numpy.abs(func(x))), i["fun"]])

                    else:
                        print("The type of constraint provided is unknown. Please review the entered constraints.")
                        sys.exit()
            else:
                self._GeneralApproach__full_constraints = []
        else:
            samples = None
            self._GeneralApproach__temp_bounds = None
            self._GeneralApproach__x_full = None
            self._GeneralApproach__full_constraints = None

        if self.__parallel_flag:
            sample_portion = mpi_mod.distribute_points(samples, self.__my_rank, self.__num_cores, self.__comm)

            # broadcast important variables
            self._GeneralApproach__temp_bounds = self.__comm.bcast(self._GeneralApproach__temp_bounds, root=0)
            self._GeneralApproach__x_full = self.__comm.bcast(self._GeneralApproach__x_full, root=0)
            self._GeneralApproach__full_constraints = self.__comm.bcast(self._GeneralApproach__full_constraints, root=0)

            self.__comm.Barrier()

            return sample_portion
        else:
            return samples

    def __initialize_optimization_mca(self):

        if self.__my_rank == 0 or self.__my_rank is None:

            # Generate starting points uniformly with length ranges
            numpy.random.seed(self._MassConservationApproach__seed)
            samples = numpy.random.rand(self.__iterations, len(self._MassConservationApproach__bounds) -
                                        len(self._MassConservationApproach__equality_bounds_indices)).astype(self._MassConservationApproach__numpy_dtype)

            x_candidates = []
            self._MassConservationApproach__x_full = numpy.zeros(len(self._MassConservationApproach__bounds), dtype=self._MassConservationApproach__numpy_dtype)

            self._MassConservationApproach__non_equality_bounds_indices = [i for i in range(len(
                self._MassConservationApproach__bounds)) if i not in self._MassConservationApproach__equality_bounds_indices]

            self._MassConservationApproach__true_bounds = [(self._MassConservationApproach__numpy_dtype(self._MassConservationApproach__bounds[j][0]),
                                                           self._MassConservationApproach__numpy_dtype(self._MassConservationApproach__bounds[j][1]))
                                                           for j in self._MassConservationApproach__non_equality_bounds_indices]

            ranges = numpy.asarray(self._MassConservationApproach__true_bounds, dtype=self._MassConservationApproach__numpy_dtype)
            samples = samples * (ranges[:, 1] - ranges[:, 0]) + ranges[:, 0]

        else:
            samples = None
            self._MassConservationApproach__x_full = None
            self._MassConservationApproach__non_equality_bounds_indices = None
            self._MassConservationApproach__true_bounds = None

        if self.__parallel_flag:
            sample_portion = mpi_mod.distribute_points(samples, self.__my_rank, self.__num_cores, self.__comm)

            # broadcast important variables
            self._MassConservationApproach__non_equality_bounds_indices = self.__comm.bcast(self._MassConservationApproach__non_equality_bounds_indices, root=0)
            self._MassConservationApproach__x_full = self.__comm.bcast(self._MassConservationApproach__x_full, root=0)
            self._MassConservationApproach__true_bounds = self.__comm.bcast(self._MassConservationApproach__true_bounds, root=0)

            self.__comm.Barrier()

            return sample_portion
        else:
            return samples

    def __initialize_optimization_sda(self):

        if self.__my_rank == 0 or self.__my_rank is None:

            # Generate starting points uniformly with length ranges
            numpy.random.seed(self._SemiDiffusiveApproach__seed)
            samples = numpy.random.rand(self.__iterations, len(self._SemiDiffusiveApproach__bounds) -
                                        len(self._SemiDiffusiveApproach__equality_bounds_indices)).astype(self._SemiDiffusiveApproach__numpy_dtype)

            x_candidates = []
            self._SemiDiffusiveApproach__x_full = numpy.zeros(len(self._SemiDiffusiveApproach__bounds), dtype=self._SemiDiffusiveApproach__numpy_dtype)

            self._SemiDiffusiveApproach__non_equality_bounds_indices = [i for i in range(len(self._SemiDiffusiveApproach__bounds))
                                                                        if i not in self._SemiDiffusiveApproach__equality_bounds_indices]

            self._SemiDiffusiveApproach__true_bounds = [(self._SemiDiffusiveApproach__numpy_dtype(self._SemiDiffusiveApproach__bounds[j][0]),
                                                        self._SemiDiffusiveApproach__numpy_dtype(self._SemiDiffusiveApproach__bounds[j][1]))
                                                        for j in self._SemiDiffusiveApproach__non_equality_bounds_indices]

            ranges = numpy.asarray(self._SemiDiffusiveApproach__true_bounds, dtype=self._SemiDiffusiveApproach__numpy_dtype)
            samples = samples * (ranges[:, 1] - ranges[:, 0]) + ranges[:, 0]

        else:
            samples = None
            self._SemiDiffusiveApproach__x_full = None
            self._SemiDiffusiveApproach__non_equality_bounds_indices = None
            self._SemiDiffusiveApproach__true_bounds = None

        if self.__parallel_flag:
            sample_portion = mpi_mod.distribute_points(samples, self.__my_rank, self.__num_cores, self.__comm)

            # broadcast important variables
            self._SemiDiffusiveApproach__non_equality_bounds_indices = self.__comm.bcast(self._SemiDiffusiveApproach__non_equality_bounds_indices, root=0)
            self._SemiDiffusiveApproach__x_full = self.__comm.bcast(self._SemiDiffusiveApproach__x_full, root=0)
            self._SemiDiffusiveApproach__true_bounds = self.__comm.bcast(self._SemiDiffusiveApproach__true_bounds, root=0)

            self.__comm.Barrier()

            return sample_portion
        else:
            return samples

    def __main_optimization_routine(self, samples):

        if self.__method == "GeneralApproach":
            feasible_point_sets = [samples[i] for i in range(len(samples))]
        else:
            feasible_point_sets = self.__feasible_point_method(samples)

        det_point_sets = []
        det_point_sets_fun = []
        smallest_value = numpy.float(1e8)
        if self.__confidence_level_flag:
            obtained_minimums = numpy.zeros(len(feasible_point_sets), dtype=numpy.float64)
        else:
            obtained_minimums = None

        if len(feasible_point_sets) != 0:

            if self.__comm is not None:
                if self.__my_rank == 0:
                    print("")
                    print("Running the multistart optimization method ...")
                    self.__start_time = MPI.Wtime()
            else:
                print("")
                print("Running the multistart optimization method ...")
                self.__start_time = time.time()

            for i in range(len(feasible_point_sets)):

                with numpy.errstate(divide='ignore', invalid='ignore'):

                    if self.__method == "GeneralApproach":
                        result = self._GeneralApproach__run_global_optimization_routine(feasible_point_sets[i])
                    elif self.__method == "MassConservationApproach":
                        result = self._MassConservationApproach__run_global_optimization_routine(feasible_point_sets[i])
                    elif self.__method == "SemiDiffusiveApproach":
                        result = self._SemiDiffusiveApproach__run_global_optimization_routine(feasible_point_sets[i])

                    if self.__print_flag:
                        print("Global function value: " + str(result.fun))
                        print("Decision vector used: ")
                        print(result.x)

                    if abs(result.fun) > numpy.float64(1e-100):

                        if self.__method == "GeneralApproach":
                            result1 = self._GeneralApproach__run_local_optimization_routine(result.x)
                        elif self.__method == "MassConservationApproach":
                            result1 = self._MassConservationApproach__run_local_optimization_routine(result.x)
                        elif self.__method == "SemiDiffusiveApproach":
                            result1 = self._SemiDiffusiveApproach__run_local_optimization_routine(result.x)

                        if self.__print_flag:
                            print("Local function value: " + str(result1.fun))
                            print("Decision vector used: ")
                            print(result1.x)

                        if smallest_value > result1.fun:
                            smallest_value = result1.fun

                        if self.__confidence_level_flag:
                            obtained_minimums[i] = result1.fun

                        if abs(result1.fun) <= numpy.finfo(float).eps:

                            if self.__method == "GeneralApproach":
                                out = self._GeneralApproach__create_final_points(result1.x)
                            elif self.__method == "MassConservationApproach":
                                out = self._MassConservationApproach__create_final_points(result1.x)
                            elif self.__method == "SemiDiffusiveApproach":
                                out = self._SemiDiffusiveApproach__create_final_points(result1.x)

                            if out is not None:
                                det_point_sets.append(out)
                                det_point_sets_fun.append(result1.fun)

                    else:
                        if smallest_value > result.fun:
                            smallest_value = result.fun

                        if self.__method == "GeneralApproach":
                            out = self._GeneralApproach__create_final_points(result.x)
                        elif self.__method == "MassConservationApproach":
                            out = self._MassConservationApproach__create_final_points(result.x)
                        elif self.__method == "SemiDiffusiveApproach":
                            out = self._SemiDiffusiveApproach__create_final_points(result.x)

                        if out is not None:
                            det_point_sets.append(out)
                            det_point_sets_fun.append(result.fun)
                        if self.__confidence_level_flag:
                            obtained_minimums[i] = result.fun

                    if self.__print_flag:
                        print("")

            if self.__comm is not None:
                self.__comm.Barrier()
                if self.__my_rank == 0:
                    self.__end_time = MPI.Wtime()
                    print("Elapsed time for multistart method: " + str(self.__end_time - self.__start_time))
                    print("")
            else:
                self.__end_time = time.time()
                print("Elapsed time for multistart method: " + str(self.__end_time - self.__start_time))
                print("")
        else:
            raise Exception("Optimization needs to be run with more iterations or different bounds.")

        return det_point_sets, det_point_sets_fun, obtained_minimums, smallest_value

    def __feasible_point_method(self, samples):

        x_candidates = []
        if self.__comm is not None:
            if self.__my_rank == 0:
                print("")
                print("Running feasible point method for " + str(self.__iterations) + " iterations ...")
                self.__start_time = MPI.Wtime()
        else:
            print("")
            print("Running feasible point method for " + str(self.__iterations) + " iterations ...")
            self.__start_time = time.time()

        for n in range(len(samples)):
            with numpy.errstate(divide='ignore', invalid='ignore'):

                if self.__method == "MassConservationApproach":
                    result = self._MassConservationApproach__run_local_optimization_routine_penalty_1(samples[n])
                elif self.__method == "SemiDiffusiveApproach":
                    result = self._SemiDiffusiveApproach__run_local_optimization_routine_penalty_1(samples[n])

                if abs(result.fun) > numpy.float64(1e-100):

                    if self.__method == "MassConservationApproach":
                        result0 = self._MassConservationApproach__run_local_optimization_routine_penalty_2(result.x)
                        output = self._MassConservationApproach__feasible_point_check(result0.x, result0.fun)
                    elif self.__method == "SemiDiffusiveApproach":
                        result0 = self._SemiDiffusiveApproach__run_local_optimization_routine_penalty_2(result.x)
                        output = self._SemiDiffusiveApproach__feasible_point_check(result0.x, result0.fun)

                    if self.__print_flag:
                        print("Objective function value: " + str(result0.fun))
                        print("Decision vector used: ")
                        print(result0.x)
                        print("")

                    if output or self.__confidence_level_flag:
                        x_candidates.append(result0.x)
                else:

                    if self.__method == "MassConservationApproach":
                        output = self._MassConservationApproach__feasible_point_check(result.x, result.fun)
                    elif self.__method == "SemiDiffusiveApproach":
                        output = self._SemiDiffusiveApproach__feasible_point_check(result.x, result.fun)

                    if self.__print_flag:
                        print("Objective function value: " + str(result.fun))
                        print("Decision vector used: ")
                        print(result.x)
                        print("")
                    if output or self.__confidence_level_flag:
                        x_candidates.append(result.x)

        if self.__comm is not None:
            self.__comm.Barrier()
            if self.__my_rank == 0:
                self.__end_time = MPI.Wtime()
                print("Elapsed time for feasible point method: " + str(self.__end_time - self.__start_time))

            # checking number of elements of feasible_point_sets for each core to see if we need to redistribute them
            redistribute_flag = len(x_candidates) == len(samples)
            val = self.__comm.allreduce(redistribute_flag, op=MPI.LAND)

            if not val:
                len_total = self.__comm.allreduce(len(x_candidates), op=MPI.SUM)
                if len_total > 0:
                    array_of_feasibles = mpi_mod.gather_numpy_array_of_values(x_candidates, self.__comm, self.__my_rank)
                    x_candidates = mpi_mod.distribute_points(array_of_feasibles, self.__my_rank, self.__num_cores,
                                                             self.__comm)
                else:
                    if self.__my_rank == 0:
                        print("No feasible points were found, please rerun the optimization routine with more iterations.")
                    sys.exit()
            self.__comm.Barrier()
        else:
            self.__end_time = time.time()
            print("Elapsed time for feasible point method: " + str(self.__end_time - self.__start_time))

            if len(x_candidates) == 0:
                print("No feasible points were found, please rerun optimization with more iterations.")
                sys.exit()

        return x_candidates

    @staticmethod
    def __confidence_level(obtained_minimums, change_in_rel_error, important_info):

        a = 1
        b = 5

        unique_elements, counts_elements = numpy.unique(obtained_minimums, return_counts=True)
        min_val_index = numpy.nanargmin(unique_elements)

        f_til = unique_elements[min_val_index]

        numpy_float64_smalles_positive_value = numpy.nextafter(numpy.float64(0), numpy.float64(1))

        if f_til > numpy_float64_smalles_positive_value:

            r = numpy.count_nonzero(
                numpy.abs(f_til - obtained_minimums) / f_til < numpy.float64(change_in_rel_error))

            n_til = obtained_minimums.shape[0]
            a_bar = a + b - 1
            b_bar = b - r - 1

            prob = 1 - (math.factorial(n_til + a_bar) * math.factorial(2 * n_til + b_bar)) / (
                    math.factorial(2 * n_til + a_bar) * math.factorial(n_til + b_bar))

        else:

            prob = 1.0

        important_info += f"It was found that {unique_elements[min_val_index]} is the minimum objective function value with a confidence level of {prob}.\n"
        return important_info
