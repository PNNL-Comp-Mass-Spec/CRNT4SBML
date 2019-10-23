from sys import platform as sys_pf
import os
import multiprocessing
import scipy
import sympy
import numpy.linalg
import time
import math
import shutil
import subprocess
import pickle
import antimony
import roadrunner
import rrplugins
import sys
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
else:
    import matplotlib.pyplot as plt


class BistabilityFinder:

    @classmethod
    def run_optimization(cls, bounds, num_constraint_method_iters, sys_min_val, temp_c, penalty_objective_func,
                         feasible_point_check, objective_function_to_optimize, final_constraint_check, seed,
                         equality_bounds_indices, print_flag, numpy_dtype, concentration_bounds):
        important_info = ''

        # running penalty method to find feasible points
        print("")
        print("Running feasible point method for " + str(num_constraint_method_iters) + " iterations ...")
        start = time.process_time()
        x_candidates = cls.feasible_point_method(bounds, num_constraint_method_iters, sys_min_val, temp_c,
                                                 penalty_objective_func, feasible_point_check, seed,
                                                 equality_bounds_indices, print_flag,
                                                 numpy_dtype, concentration_bounds)
        end = time.process_time()
        print("Elapsed time for feasible point method: " + str(end - start))
        print("")

        num_x_cand = len(x_candidates)
        important_info += "The number of feasible points that satisfy the constraints: " + str(num_x_cand) + "\n"

        # if num_x_cand = 0 throw error and don't cont.
        if num_x_cand != 0:
            print("Running the multistart optimization ...")
            start = time.process_time()

            params_that_give_global_min, total_that_give_zero, total_that_pass_chk, obj_fun_val_global_min = \
                cls.multistart_optimization(num_x_cand, x_candidates, sys_min_val, bounds, temp_c,
                                            objective_function_to_optimize, equality_bounds_indices, seed,
                                            final_constraint_check, print_flag, numpy_dtype, concentration_bounds)
            
            end = time.process_time()
            print("Elapsed time for multistart method: " + str(end - start))
            print("")

            important_info += "Total feasible points that give F(x) = 0: " + str(total_that_give_zero) + "\n"

            important_info += "Total number of points that passed final_check: " + str(total_that_pass_chk) + "\n"

            return params_that_give_global_min, obj_fun_val_global_min, important_info
        else:
            raise Exception("Optimization needs to be run with more iterations or different bounds.")

    @classmethod
    def feasible_point_method(cls, penalty_bounds, num_constraint_method_iters, sys_min_val, temp_c,
                              penalty_objective_func, feasible_point_check, seed, equality_bounds_indices, print_flag,
                              numpy_dtype, concentration_bounds):
        # Generate starting points uniformly with length ranges
        numpy.random.seed(seed)
        samples = numpy.random.rand(num_constraint_method_iters, len(penalty_bounds) -
                                    len(equality_bounds_indices)).astype(numpy_dtype)

        x_candidates = []
        x_off = numpy.zeros(len(penalty_bounds), dtype=numpy_dtype)

        non_equality_bounds_indices = [i for i in range(len(penalty_bounds)) if i not in equality_bounds_indices]
        
        true_bounds = [(numpy_dtype(penalty_bounds[j][0]), numpy_dtype(penalty_bounds[j][1]))
                       for j in non_equality_bounds_indices]

        ranges = numpy.asarray(true_bounds, dtype=numpy_dtype)
        samples = samples*(ranges[:, 1] - ranges[:, 0]) + ranges[:, 0]

        for n in range(num_constraint_method_iters):
            with numpy.errstate(divide='ignore', invalid='ignore'):
                result = scipy.optimize.minimize(penalty_objective_func, samples[n], args=
                (temp_c, penalty_bounds, equality_bounds_indices, x_off, non_equality_bounds_indices,
                 concentration_bounds), method='SLSQP', tol=1e-16, bounds=true_bounds)

                if abs(result.fun) > numpy_dtype(1e-100): 
                    result0 = scipy.optimize.minimize(penalty_objective_func, result.x, args=
                    (temp_c, penalty_bounds, equality_bounds_indices, x_off, non_equality_bounds_indices,
                     concentration_bounds), method='Nelder-Mead', tol=1e-16)

                    output = feasible_point_check(result0.x, result0.fun, sys_min_val, equality_bounds_indices,
                                                  non_equality_bounds_indices, penalty_bounds, concentration_bounds)
                    if print_flag:
                        print("Objective function value: " + str(result0.fun))
                        print("Decision vector used: ")
                        print(result0.x)
                        print("")

                    if output:
                        x_candidates.append(result0.x)
                else:
                    output = feasible_point_check(result.x, result.fun, sys_min_val, equality_bounds_indices,
                                                  non_equality_bounds_indices, penalty_bounds, concentration_bounds)
                    if print_flag:
                        print("Objective function value: " + str(result.fun))
                        print("Decision vector used: ")
                        print(result.x)
                        print("")
                    if output:
                        x_candidates.append(result.x)

        return x_candidates

    @classmethod 
    def multistart_optimization(cls, num_x_cand, x_candidates, sys_min_val, penalty_bounds, temp_c,
                                objective_function_to_optimize, equality_bounds_indices, seed_given,
                                final_constraint_check, print_flag, numpy_dtype, concentration_bounds):
        x_off = numpy.zeros(len(penalty_bounds), dtype=numpy_dtype)
        non_equality_bounds_indices = [i for i in range(len(penalty_bounds)) if i not in equality_bounds_indices]

        smallest_value = numpy_dtype(numpy.PINF)
        x_that_give_global_min = []
        obj_fun_val_global_min = []
        start = 0
        for i in range(num_x_cand):
            with numpy.errstate(divide='ignore', invalid='ignore'):
                result = scipy.optimize.basinhopping(objective_function_to_optimize, x_candidates[i],
                                                     minimizer_kwargs={'method': 'Nelder-Mead',
                                                                       'args': (temp_c, penalty_bounds, sys_min_val,
                                                                                equality_bounds_indices, x_off,
                                                                                non_equality_bounds_indices,
                                                                                concentration_bounds), 'tol': 1e-16},
                                                     niter=2, seed=seed_given)

                start = start+1

                if print_flag:
                    print("Global function value: " + str(result.fun))
                    print("Decision vector produced: ")
                    print(result.x)
                
                if result.fun > numpy_dtype(1e-200):
                    result = scipy.optimize.minimize(objective_function_to_optimize, result.x,
                                                     args=(temp_c, penalty_bounds, sys_min_val, equality_bounds_indices,
                                                           x_off, non_equality_bounds_indices, concentration_bounds),
                                                     method='Nelder-Mead', tol=1e-16)
        
                    if print_flag:
                        print("Local function value: " + str(result.fun))
                        print("Decision vector produced: ")
                        print(result.x)

                if print_flag:
                    print("")
                        
                if result.fun < smallest_value:
                    smallest_value = result.fun

            if abs(result.fun) <= sys_min_val: 
                x_that_give_global_min.append(result.x)
                obj_fun_val_global_min.append(result.fun)

        x_that_give_global_min2, total_that_give_zero, total_that_pass_chk, obj_fun_val_global_min2 = \
            cls.__create_final_points(x_that_give_global_min, obj_fun_val_global_min, final_constraint_check,
                                      penalty_bounds, sys_min_val, equality_bounds_indices, concentration_bounds)

        print("\nSmallest value achieved by objective function: " + str(smallest_value) + "\n")

        return x_that_give_global_min2, total_that_give_zero, total_that_pass_chk, obj_fun_val_global_min2

    @staticmethod
    def __create_final_points(x_that_give_global_min, obj_fun_val_global_min, final_constraint_check, penalty_bounds,
                              sys_min_val, non_equality_bounds_indices, concentration_bounds):
        if len(x_that_give_global_min) != 0:
            temp2 = numpy.vstack(x_that_give_global_min)
            u, indices = numpy.unique(temp2, axis=0, return_index=True)
            x_that_give_global_min2 = []
            obj_fun_val_global_min2 = []
            for i in indices:

                output = final_constraint_check(x_that_give_global_min[i], penalty_bounds, sys_min_val,
                                                non_equality_bounds_indices, concentration_bounds)

                if output[0]:
                    x_that_give_global_min2.append(output[1])
                    obj_fun_val_global_min2.append(obj_fun_val_global_min[i])
            total_that_give_zero = len(indices)
            total_that_pass_chk = len(x_that_give_global_min2)
        else:
            x_that_give_global_min2 = []
            obj_fun_val_global_min2 = []
            total_that_give_zero = 0
            total_that_pass_chk = 0

        return x_that_give_global_min2, total_that_give_zero, total_that_pass_chk, obj_fun_val_global_min2
        
    @classmethod
    def run_continuity_analysis(cls, species_num, params_for_global_min, initialize_ant_string, finalize_ant_string,
                                species_y, dir_path, print_lbls_flag, auto_parameters):

        print("Running continuity analysis ...")

        if sys_pf not in ['win32', 'cygwin', 'msys']:
            roadrunner.Logger.setLevel(roadrunner.Logger.LOG_ERROR)
            roadrunner.Logger.disableLogging()
            roadrunner.Logger.disableConsoleLogging()
            roadrunner.Logger.disableFileLogging()
            rrplugins.setLogLevel('error')
            try:
                stderr_fileno = sys.stderr.fileno()
                stderr_save = os.dup(stderr_fileno)
                stderr_pipe = os.pipe()
                os.dup2(stderr_pipe[1], stderr_fileno)
                os.close(stderr_pipe[1])
                notebook_exists = False
            except Exception as e:
                print("Note: stderr is not being caught in the traditional fashion. This may be a result of using a notebook.")
                notebook_exists = True

        init_ant, pcp_x = initialize_ant_string(species_num, auto_parameters['PrincipalContinuationParameter'])
        auto_parameters['PrincipalContinuationParameter'] = pcp_x
        start = time.time()
        multistable_param_ind = []
        cont_direction = ["Positive", "Negative"]

        if sys_pf not in ['win32', 'cygwin', 'msys']:
            auto = rrplugins.Plugin("tel_auto2000")
        else:
            auto = None

        plot_specifications = []

        for param_ind in range(len(params_for_global_min)):
            final_ant_str = finalize_ant_string(params_for_global_min[param_ind], init_ant)

            for dir_ind in range(2):
                if os.path.isdir("./auto_fort_files"):
                    shutil.rmtree("./auto_fort_files")

                pts, lbls, antimony_r, flag, bi_data_np = cls.run_safety_wrapper(final_ant_str, cont_direction[dir_ind],
                                                                                 auto, auto_parameters)

                if print_lbls_flag:
                    print("Labels from numerical continuation: ")
                    print(lbls)
                
                if flag and lbls != ['EP', 'EP']:
                    chnk_stable, chnk_unstable, special_points, sp_y_ind = cls.find_stable_unstable_regions(antimony_r,
                                                                                                            species_y)

                    multistable = cls.detect_multi_stability(chnk_stable, chnk_unstable, bi_data_np)

                    if multistable:
                        plt_specs = cls.plot_pcp_vs_species(chnk_stable, chnk_unstable, special_points, bi_data_np,
                                                            sp_y_ind, pcp_x, species_y, param_ind, dir_path, cls)
                        plot_specifications.append(plt_specs)
                        multistable_param_ind.append(param_ind)
                        break

        if os.path.isdir("./auto_fort_files"):
            shutil.rmtree("./auto_fort_files")

        if sys_pf not in ['win32', 'cygwin', 'msys']:
            if not notebook_exists:
                os.close(stderr_pipe[0])
                os.dup2(stderr_save, stderr_fileno)
                os.close(stderr_save)
                os.close(stderr_fileno)

        end = time.time()
        print("Elapsed time for continuity analysis: " + str(end - start))
        print("")

        important_info = "Number of multistability plots found: " + str(len(multistable_param_ind)) + "\n"
        
        important_info += "Elements in params_for_global_min that produce multistability: \n" + \
                          str(multistable_param_ind) + "\n"

        return multistable_param_ind, important_info, plot_specifications

    @classmethod
    def run_greedy_continuity_analysis(cls, species_num, params_for_global_min, initialize_ant_string,
                                       finalize_ant_string,
                                       species_y, dir_path, print_lbls_flag, auto_parameters):

        print("Running continuity analysis ...")

        if sys_pf not in ['win32', 'cygwin', 'msys']:
            roadrunner.Logger.setLevel(roadrunner.Logger.LOG_ERROR)
            roadrunner.Logger.disableLogging()
            roadrunner.Logger.disableConsoleLogging()
            roadrunner.Logger.disableFileLogging()
            rrplugins.setLogLevel('error')
            try:
                sys.stderr.fileno()
                stderr_fileno = sys.stderr.fileno()
                stderr_save = os.dup(stderr_fileno)
                stderr_pipe = os.pipe()
                os.dup2(stderr_pipe[1], stderr_fileno)
                os.close(stderr_pipe[1])
                notebook_exists = False
            except Exception as e:
                print("Note: stderr is not being caught in the traditional fashion. This may be a result of using a notebook.")
                notebook_exists = True

        init_ant, pcp_x = initialize_ant_string(species_num, auto_parameters['PrincipalContinuationParameter'])
        auto_parameters['PrincipalContinuationParameter'] = pcp_x
        auto_parameters['IADS'] = 0
        auto_parameters['A1'] = 1e10
        auto_parameters['ITNW'] = 100
        auto_parameters['NTST'] = 100
        auto_parameters['NCOL'] = 100

        start = time.time()
        multistable_param_ind = []
        cont_direction = ["Positive", "Negative"]

        if sys_pf not in ['win32', 'cygwin', 'msys']:
            auto = rrplugins.Plugin("tel_auto2000")
        else:
            auto = None

        plot_specifications = []

        for param_ind in range(len(params_for_global_min)):

            final_ant_str = finalize_ant_string(params_for_global_min[param_ind], init_ant)
            pcp_ranges_mag = cls.get_pcp_ranges(final_ant_str, pcp_x)

            auto_parameters['RL0'] = pcp_ranges_mag[0]
            auto_parameters['RL1'] = pcp_ranges_mag[1]

            ds_vals = []
            mag = pcp_ranges_mag[2]
            mag -= 1
            for i in range(5):
                ds_vals.append(float(10 ** mag))
                mag -= 1

            multistable = False
            lbls = []
            for ds_val in ds_vals:
                auto_parameters['DS'] = ds_val
                for dir_ind in range(2):
                    if os.path.isdir("./auto_fort_files"):
                        shutil.rmtree("./auto_fort_files")

                    pts, lbls, antimony_r, flag, bi_data_np = cls.run_safety_wrapper(final_ant_str,
                                                                                     cont_direction[dir_ind], auto,
                                                                                     auto_parameters)

                    if print_lbls_flag:
                        print("Labels from numerical continuation: ")
                        print(lbls)

                    if flag and lbls != ['EP', 'EP']:
                        chnk_stable, chnk_unstable, special_points, sp_y_ind = cls.find_stable_unstable_regions(
                            antimony_r,
                            species_y)

                        multistable = cls.detect_multi_stability(chnk_stable, chnk_unstable, bi_data_np)

                        if multistable:

                            # running another numerical continuation with a smaller step size to try and get
                            # better looking plots
                            scnd_check = True
                            if os.path.isdir("./auto_fort_files"):
                                shutil.rmtree("./auto_fort_files")

                                ind_ds_val = ds_vals.index(ds_val)
                                ds_val = ds_vals[ind_ds_val + 1]

                                auto_parameters['DS'] = ds_val
                                for dir_ind2 in range(2):
                                    if os.path.isdir("./auto_fort_files"):
                                        shutil.rmtree("./auto_fort_files")

                                    pts2, lbls2, antimony_r, flag2, bi_data_np2 = cls.run_safety_wrapper(final_ant_str,
                                                                                                         cont_direction[dir_ind2],
                                                                                                         auto,
                                                                                                         auto_parameters)

                                    if flag2 and lbls2 != ['EP', 'EP']:
                                        chnk_stable2, chnk_unstable2, special_points2, sp_y_ind2 = \
                                            cls.find_stable_unstable_regions(antimony_r, species_y)

                                        multistable2 = cls.detect_multi_stability(chnk_stable2, chnk_unstable2,
                                                                                  bi_data_np2)

                                        if multistable2:
                                            plt_specs = cls.plot_pcp_vs_species(chnk_stable2, chnk_unstable2,
                                                                                special_points2, bi_data_np2, sp_y_ind2,
                                                                                pcp_x, species_y, param_ind, dir_path,
                                                                                cls)
                                            plot_specifications.append(plt_specs)
                                            scnd_check = False
                                            break

                                if scnd_check:
                                    plt_specs = cls.plot_pcp_vs_species(chnk_stable, chnk_unstable, special_points,
                                                                        bi_data_np, sp_y_ind, pcp_x, species_y,
                                                                        param_ind, dir_path, cls)
                                    plot_specifications.append(plt_specs)

                            if param_ind not in multistable_param_ind:
                                multistable_param_ind.append(param_ind)
                            break
                if multistable and 'MX' not in lbls:
                    break
            if print_lbls_flag:
                print("")

        if os.path.isdir("./auto_fort_files"):
            shutil.rmtree("./auto_fort_files")

        if sys_pf not in ['win32', 'cygwin', 'msys']:
            if not notebook_exists:
                os.close(stderr_pipe[0])
                os.dup2(stderr_save, stderr_fileno)
                os.close(stderr_save)
                os.close(stderr_fileno)

        end = time.time()
        print("Elapsed time for continuity analysis: " + str(end - start))
        print("")

        important_info = "Number of multistability plots found: " + str(len(multistable_param_ind)) + "\n"

        important_info += "Elements in params_for_global_min that produce multistability: \n" + \
                          str(multistable_param_ind) + "\n"

        return multistable_param_ind, important_info, plot_specifications

    @staticmethod
    def get_pcp_ranges(final_ant_str, pcp_x):

        splits = final_ant_str.split(';')
        initial_vars = [i for i in splits if '=' in i]
        symbol_and_vals = []
        for j in initial_vars:
            symbol_and_vals.append([i.strip() for i in j.split('=')])

        pcp_initial = [symbol_and_vals[i] for i in range(len(symbol_and_vals))
                       if symbol_and_vals[i][0] == pcp_x][0][1]

        def order_of_magnitude(number):
            return math.floor(math.log(number, 10))

        mag = order_of_magnitude(float(pcp_initial))

        for i in pcp_initial:
            if i != '0' and i != '.':
                val = i
                break

        rl1 = float(val) * 10 ** (mag + 1)
        rl0 = float(val) * 10 ** (mag - 1)

        return [rl0, rl1, mag + 1]

    @classmethod
    def run_safety_wrapper(cls, final_ant_str, cont_direction, auto, auto_parameters):

        if sys_pf in ['win32', 'cygwin', 'msys']:
            arguments = [final_ant_str, cont_direction, auto_parameters]
            if os.path.exists("input_arguments.pickle"):
                os.remove("input_arguments.pickle")
                with open('input_arguments.pickle', 'wb') as outf:
                    outf.write(pickle.dumps(arguments))
            else:
                with open('input_arguments.pickle', 'wb') as outf:
                    outf.write(pickle.dumps(arguments))

            # # making the directory auto_fort_files if is does not exist
            if not os.path.isdir("./auto_fort_files"):
                os.mkdir("./auto_fort_files")

            subprocess.run(['python', '-m', 'crnt4sbml.safety_wrap'], shell=True, env=os.environ)

            if os.path.exists("output_arguments.pickle"):
                with open('output_arguments.pickle', 'rb') as pickle_file:
                    output_arguments = pickle.loads(pickle_file.read())
                os.remove("output_arguments.pickle")
                pts = output_arguments[0]
                lbls = output_arguments[1]
                antimony_r = output_arguments[2]
                flag = output_arguments[3]
                bi_data_np = numpy.load('bi_data_np.npy')
                os.remove('./bi_data_np.npy')
                if os.path.exists("input_arguments.pickle"):
                    os.remove("input_arguments.pickle")
            else:
                flag = False
                pts = []
                lbls = []
                antimony_r = []
                bi_data_np = numpy.zeros(1)
                if os.path.exists("input_arguments.pickle"):
                    os.remove("input_arguments.pickle")

            time.sleep(1)
        else:

            if os.path.exists("bi_data_np.npy"):
                os.remove("bi_data_np.npy")

            queue = multiprocessing.Queue()
            p = multiprocessing.Process(target=cls.run_numerical_continuation, args=(queue, final_ant_str,
                                                                                     cont_direction,
                                                                                     auto, auto_parameters))

            p.start()
            p.join()  # this blocks until the process terminates

            if p.exitcode == 0:
                pts, lbls, antimony_r, flag = queue.get()
                bi_data_np = numpy.load('bi_data_np.npy')
                if os.path.exists("bi_data_np.npy"):
                    os.remove("bi_data_np.npy")
                p.terminate()
                queue.close()
            else:
                flag = False
                pts = []
                lbls = []
                antimony_r = []
                bi_data_np = numpy.zeros(1)
                p.terminate()
                queue.close()
                if os.path.exists("bi_data_np.npy"):
                    os.remove("bi_data_np.npy")

        return pts, lbls, antimony_r, flag, bi_data_np

    @classmethod
    def run_numerical_continuation(cls, q, ant_str, direction, auto, auto_parameters):

        antimony_r = cls.__loada(ant_str)

        if True:

            # making the directory auto_fort_files if is does not exist
            if not os.path.isdir("./auto_fort_files"):
                os.mkdir("./auto_fort_files")

            auto.setProperty("SBML", antimony_r.getCurrentSBML())
            auto.setProperty("ScanDirection", direction)
            auto.setProperty("PreSimulation", "True")
            auto.setProperty("PreSimulationDuration", 1.0)
            auto.setProperty('KeepTempFiles', True)
            auto.setProperty("TempFolder", "./auto_fort_files")

            # assigning values provided by the user
            for i in auto_parameters.keys():
                auto.setProperty(i, auto_parameters[i])

            try:
                auto.execute()
                # indices where special points are
                pts = auto.BifurcationPoints
                # labeling of special points
                lbls = auto.BifurcationLabels
                # all data for parameters and species found by continuation
                bi_data = auto.BifurcationData

                # convertes bi_data to numpy array, where first
                # column is the principal continuation parameter and
                # the rest of the columns are the species
                bi_data_np = bi_data.toNumpy
                flag = True

            except Exception as e:
                flag = False
                pts = []
                lbls = []
                bi_data_np = numpy.zeros(2)

        ant_float_ids = antimony_r.model.getFloatingSpeciesIds()
        numpy.save('bi_data_np.npy', bi_data_np)
        q.put([pts, lbls, ant_float_ids, flag])

    @staticmethod
    def find_stable_unstable_regions(antimony_r, species_y):
        with open("./auto_fort_files/fort.7", 'r') as fobj:
            all_lines = [[num for num in line.split()] for line in fobj]
        
        # getting index where important information begins
        beg_ind = None
        for i in range(len(all_lines)):
            if all_lines[i][0] == '1':
                beg_ind = i
                break

        solution_types = {'1': 'BP', '2': 'LP', '3': 'HB', '5': 'LP', '6': 'BP', '7': 'PD', '8': 'TR',
                          '9': 'EP', '-9': 'MX'}

        # lists to hold stable and unstable regions, and special points
        unstable = []
        stable = []
        special_points = []
        ignored_labels = ['0', '4', '-4']
        for i in range(beg_ind, len(all_lines)):
            if all_lines[i][0] == '0':
                break
            else:
                temp = all_lines[i][1:3]
                if temp[1] in ignored_labels:
                    if temp[0][0] == "-":
                        stable.append(abs(int(temp[0]))-1)
                    else:
                        unstable.append(int(temp[0])-1)
                else:
                    special_points.append((abs(int(temp[0]))-1, solution_types[temp[1]]))
            
        # getting stable chunks of indices
        spl = [0] + [i for i in range(1, len(stable)) if stable[i] - stable[i-1] > 1]+ [None]
        chnk_stable = [stable[b:e] for (b, e) in [(spl[i-1], spl[i]) for i in range(1, len(spl))]]

        # getting unstable chunks of indices
        spl = [0] + [i for i in range(1, len(unstable)) if unstable[i] - unstable[i-1] > 1] + [None]
        chnk_unstable = [unstable[b:e] for (b, e) in [(spl[i-1], spl[i]) for i in range(1, len(spl))]]

        # species ids according to AUTO also order in biData
        spec_id = antimony_r

        # index of species_y in biData
        sp_y_ind = spec_id.index(species_y) + 1

        return chnk_stable, chnk_unstable, special_points, sp_y_ind

    @classmethod
    def detect_multi_stability(cls, chnk_stable, chnk_unstable, bi_data_np):

        if cls.is_list_empty(chnk_stable) or cls.is_list_empty(chnk_unstable):
            return False

        chnk_stable_pcp_ranges = []
        for i in range(len(chnk_stable)):
            end_ind = len(chnk_stable[i])-1 
            chnk_stable_pcp_ranges.append([bi_data_np[chnk_stable[i][0], 0], bi_data_np[chnk_stable[i][end_ind], 0]])

        [i.sort() for i in chnk_stable_pcp_ranges] 

        chnk_unstable_pcp_ranges = []
        for i in range(len(chnk_unstable)):
            end_ind = len(chnk_unstable[i])-1
            chnk_unstable_pcp_ranges.append([bi_data_np[chnk_unstable[i][0], 0],
                                             bi_data_np[chnk_unstable[i][end_ind], 0]])

        [i.sort() for i in chnk_unstable_pcp_ranges] 

        chnk_unstable_pcp_intervals = []
        chnk_stable_pcp_intervals = []
        for i in range(len(chnk_unstable)):
            chnk_unstable_pcp_intervals.append(sympy.Interval(chnk_unstable_pcp_ranges[i][0],
                                                              chnk_unstable_pcp_ranges[i][1]))

        for i in range(len(chnk_stable)):
            chnk_stable_pcp_intervals.append(sympy.Interval(chnk_stable_pcp_ranges[i][0],
                                                            chnk_stable_pcp_ranges[i][1]))

        # building intersections of unstable branch with stable branches
        unstable_intersections = []
        for i in chnk_unstable_pcp_intervals:
            temp_list = []
            for j in chnk_stable_pcp_intervals:
                temp = i.intersect(j)
                if temp != sympy.EmptySet():

                    if not temp.is_FiniteSet:
                        temp_list.append(1)
                    elif temp.is_FiniteSet and len(list(temp)) > 1:
                        temp_list.append(1)

            unstable_intersections.append(temp_list)

        return any([sum(i) >= 2 for i in unstable_intersections])

    @classmethod
    def is_list_empty(cls, in_list):
        if isinstance(in_list, list):  # Is a list
            return all(map(cls.is_list_empty, in_list))
        return False  # Not a list

    @staticmethod
    def plot_pcp_vs_species(chnk_stable, chnk_unstable, special_points, bi_data_np, sp_y_ind, pcp_x, species_y,
                            param_ind, dir_path, cls):

        # plotting stable points
        for i in range(len(chnk_stable)):
            plt.plot(bi_data_np[chnk_stable[i], 0], bi_data_np[chnk_stable[i], sp_y_ind], linewidth=1,
                     linestyle='-', color='b')

        # plotting unstable points
        for i in range(len(chnk_unstable)):
            plt.plot(bi_data_np[chnk_unstable[i], 0], bi_data_np[chnk_unstable[i], sp_y_ind], linewidth=1,
                     linestyle='--', color='b')

        # plotting the special points with red marker *
        for i in special_points:
            plt.plot(bi_data_np[i[0], 0], bi_data_np[i[0], sp_y_ind], marker='*', color='r')
            plt.annotate(i[1], (bi_data_np[i[0], 0], bi_data_np[i[0], sp_y_ind]))

        margin = 1e-6
        plt.xlabel(pcp_x)
        plt.ylabel(species_y)

        pcp_values = []
        for i in special_points:
            if i[1] != 'EP':
                pcp_values.append(bi_data_np[i[0], 0])

        pcp_min = min(pcp_values)
        pcp_max = max(pcp_values)

        if abs(pcp_min - pcp_max) > 1e-16:
            diff_x = pcp_max - pcp_min
            start = pcp_min - diff_x
            end = pcp_max + diff_x
            xlim_list = [start, end]
            plt.xlim(xlim_list)
        else:
            start = pcp_min - pcp_min*margin
            end = pcp_max + pcp_max*margin
            xlim_list = [start, end]
            plt.xlim(xlim_list)

        lims = plt.gca().get_xlim()
        x = bi_data_np[:, 0]
        y = bi_data_np[:, sp_y_ind]
        i = numpy.where((x > lims[0]) & (x < lims[1]))[0]

        h = y[i].max() - y[i].min()
        first = y[i].min() - 0.1*h
        second = y[i].max() + 0.1*h
        plt.gca().set_ylim(first, second)

        plt.ticklabel_format(axis='both', style='sci', scilimits=(-2, 2))
        plt.savefig(dir_path + '/' + pcp_x + '_vs_' + species_y + '_' + str(param_ind) + '.png')
        plt.clf()

        plot_specs = [[lims[0], lims[1]], [first, second]]
        special_pnts_in_plot = []
        for iii in special_points:
            xval = bi_data_np[iii[0], 0]
            yval = bi_data_np[iii[0], sp_y_ind]
            if (lims[0] <= xval) and (lims[1] >= xval) and (first <= yval) and (second >= yval):
                special_pnts_in_plot.append([xval, yval, iii[1]])

        plot_specs.append(special_pnts_in_plot)

        return plot_specs

    # functions taken from Tellurium!! Give them
    # credit, they deserve it!
    #################################################
    @staticmethod
    def __check_antimony_return_code(code):

        if code < 0:
            raise Exception('Antimony: {}'.format(antimony.getLastError()))

    @classmethod
    def __antimony_to_sbml(cls, ant):
        try:
            isfile = os.path.isfile(ant)
        except ValueError:
            isfile = False
        if isfile:
            code = antimony.loadAntimonyFile(ant)
        else:
            code = antimony.loadAntimonyString(ant)
        cls.__check_antimony_return_code(code)
        mid = antimony.getMainModuleName()
        return antimony.getSBMLString(mid)

    @classmethod
    def __loada(cls, ant):

        return cls.__load_antimony_model(ant)

    @classmethod
    def __load_antimony_model(cls, ant):

        sbml = cls.__antimony_to_sbml(ant)
        return roadrunner.RoadRunner(sbml)
