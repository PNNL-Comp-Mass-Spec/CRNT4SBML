
class BistabilityAnalysis:

    def __init__(self, GA):

        self.__GA = GA

        print(GA)
        print(BistabilityAnalysis.__mro__)
        print(GA.__)
        print(self.__GA.__sympy_species)

    @classmethod
    def run_continuity_analysis(cls, species_num, params_for_global_min, initialize_ant_string, finalize_ant_string,
                                species_y, dir_path, print_lbls_flag, auto_parameters, plot_labels):

        global sys_pf
        global os
        global multiprocessing
        global numpy
        global sympy
        global time
        global math
        global shutil
        global subprocess
        global pickle
        global antimony
        global roadrunner
        global rrplugins
        global sys
        global plt

        from sys import platform as sys_pf
        import os
        import multiprocessing
        import sympy
        import numpy
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
                print(
                    "Note: stderr is not being caught in the traditional fashion. This may be a result of using a notebook.")
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
                                                            sp_y_ind, pcp_x, species_y, param_ind, dir_path, cls,
                                                            plot_labels)
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
        print("Elapsed time for continuity analysis in seconds: " + str(end - start))
        print("")

        important_info = "Number of multistability plots found: " + str(len(multistable_param_ind)) + "\n"

        important_info += "Elements in params_for_global_min that produce multistability: \n" + \
                          str(multistable_param_ind) + "\n"

        return multistable_param_ind, important_info, plot_specifications

    @classmethod
    def run_greedy_continuity_analysis(cls, species_num, params_for_global_min, initialize_ant_string,
                                       finalize_ant_string, species_y, dir_path, print_lbls_flag, auto_parameters, plot_labels):

        global sys_pf
        global os
        global multiprocessing
        global sympy
        global numpy
        global time
        global math
        global shutil
        global subprocess
        global pickle
        global antimony
        global roadrunner
        global rrplugins
        global sys
        global plt

        from sys import platform as sys_pf
        import os
        import multiprocessing
        import sympy
        import numpy
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
                print(
                    "Note: stderr is not being caught in the traditional fashion. This may be a result of using a notebook.")
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
                        # print(multistable)
                        # print(ds_val)
                        # print(dir_ind)
                        if multistable:  ################################################
                            # if ds_val == 0.01 and dir_ind == 0:
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
                                                                                                         cont_direction[
                                                                                                             dir_ind2],
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
                                                                                cls, plot_labels)
                                            plot_specifications.append(plt_specs)
                                            scnd_check = False
                                            break

                                if scnd_check:
                                    plt_specs = cls.plot_pcp_vs_species(chnk_stable, chnk_unstable, special_points,
                                                                        bi_data_np, sp_y_ind, pcp_x, species_y,
                                                                        param_ind, dir_path, cls, plot_labels)
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
        print("Elapsed time for continuity analysis in seconds: " + str(end - start))
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

        mag = order_of_magnitude(abs(float(pcp_initial)))

        for i in pcp_initial:
            if i != '0' and i != '.' and i!= '-':
                val = i
                break

        if pcp_initial[0] == '-':
            rl0 = -float(val) * 10 ** (mag + 1)
            rl1 = -float(val) * 10 ** (mag - 1)
            #rl0 = -100.0
            rl1 = 10.0
        else:
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
    def run_numerical_continuation(cls, q, ant_str, direction, auto, auto_parameters, core=None):

        antimony_r = cls.__loada(ant_str)

        # making the directory auto_fort_files if is does not exist
        if core is None:
            if not os.path.isdir("./auto_fort_files"):
                os.mkdir("./auto_fort_files")
        else:
            if not os.path.isdir("./auto_fort_files_" + str(core)):
                os.mkdir("./auto_fort_files_" + str(core))

        auto.setProperty("SBML", antimony_r.getCurrentSBML())
        auto.setProperty("ScanDirection", direction)
        auto.setProperty("PreSimulation", "True")
        auto.setProperty("PreSimulationDuration", 1.0)
        auto.setProperty('KeepTempFiles', True)
        if core is None:
            auto.setProperty("TempFolder", "./auto_fort_files")
        else:
            auto.setProperty("TempFolder", "./auto_fort_files_" + str(core))

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
        if core is None:
            numpy.save('bi_data_np.npy', bi_data_np)
        else:
            numpy.save('bi_data_np_' + str(core) + '.npy', bi_data_np)
        q.put([pts, lbls, ant_float_ids, flag])

    @staticmethod
    def find_stable_unstable_regions(antimony_r, species_y, core=None):

        if core is None:
            with open("./auto_fort_files/fort.7", 'r') as fobj:
                all_lines = [[num for num in line.split()] for line in fobj]
        else:
            with open("./auto_fort_files" + "_" + str(core) + "/fort.7", 'r') as fobj:
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
                        stable.append(abs(int(temp[0])) - 1)
                    else:
                        unstable.append(int(temp[0]) - 1)
                else:
                    special_points.append((abs(int(temp[0])) - 1, solution_types[temp[1]]))

        # getting stable chunks of indices
        spl = [0] + [i for i in range(1, len(stable)) if stable[i] - stable[i - 1] > 1] + [None]
        chnk_stable = [stable[b:e] for (b, e) in [(spl[i - 1], spl[i]) for i in range(1, len(spl))]]

        # getting unstable chunks of indices
        spl = [0] + [i for i in range(1, len(unstable)) if unstable[i] - unstable[i - 1] > 1] + [None]
        chnk_unstable = [unstable[b:e] for (b, e) in [(spl[i - 1], spl[i]) for i in range(1, len(spl))]]

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
            end_ind = len(chnk_stable[i]) - 1
            chnk_stable_pcp_ranges.append([bi_data_np[chnk_stable[i][0], 0], bi_data_np[chnk_stable[i][end_ind], 0]])

        [i.sort() for i in chnk_stable_pcp_ranges]

        chnk_unstable_pcp_ranges = []
        for i in range(len(chnk_unstable)):
            end_ind = len(chnk_unstable[i]) - 1
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
                            param_ind, dir_path, cls, plot_labels, core=None):

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
        if plot_labels != None:

            if plot_labels[0] != None:
                plt.xlabel(plot_labels[0])
            else:
                plt.xlabel(pcp_x)

            if plot_labels[1] != None:
                plt.ylabel(plot_labels[1])
            else:
                plt.ylabel(species_y)

            if plot_labels[2] != None:
                plt.title(plot_labels[2])

        else:
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
            start = pcp_min - pcp_min * margin
            end = pcp_max + pcp_max * margin
            xlim_list = [start, end]
            plt.xlim(xlim_list)

        lims = plt.gca().get_xlim()
        x = bi_data_np[:, 0]
        y = bi_data_np[:, sp_y_ind]
        i = numpy.where((x > lims[0]) & (x < lims[1]))[0]

        h = y[i].max() - y[i].min()
        first = y[i].min() - 0.1 * h
        second = y[i].max() + 0.1 * h
        plt.gca().set_ylim(first, second)

        plt.ticklabel_format(axis='both', style='sci', scilimits=(-2, 2))

        if core is None:
            plt.savefig(dir_path + '/' + pcp_x + '_vs_' + species_y + '_' + str(param_ind) + '.png')
            plt.clf()
        else:
            plt.savefig(
                dir_path + '/' + pcp_x + '_vs_' + species_y + '_' + str(param_ind) + "_core_" + str(core) + '.png')
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

    @classmethod
    def __print_initial_conditions(cls, fwd_scan_vals, rvrs_scan_vals):

        fwd_spec_inds = [i[0] for i in fwd_scan_vals]
        init_vals = []
        for i in range(cls.__N):
            if i in fwd_spec_inds:
                init_vals.append(str(cls.__sympy_species[i]) + " = " + "C" + str(fwd_spec_inds.index(i) + 1))
            else:
                init_vals.append(str(cls.__sympy_species[i]) + " = 0.0")

        print(" ")
        print("For the forward scan the following initial condition will be used:")
        for i in init_vals:
            print(i)

        rvrs_spec_inds = [i[0] for i in rvrs_scan_vals]
        init_vals = []
        for i in range(cls.__N):
            if i in rvrs_spec_inds:
                init_vals.append(str(cls.__sympy_species[i]) + " = " + "C" + str(rvrs_spec_inds.index(i) + 1))
            else:
                init_vals.append(str(cls.__sympy_species[i]) + " = 0.0")

        print(" ")
        print("For the reverse scan the following initial condition will be used:")
        for i in init_vals:
            print(i)
        print(" ")

    @classmethod
    def construct_conservation_vals(cls, params_for_global_min):

        if cls.__approach == "GeneralApproach":

            return [cls.__cons_laws_sympy_lamb[ii](*tuple(params_for_global_min[cls.__R:cls.__R + cls.__N]))
                    for ii in range(len(cls.__cons_laws_sympy_lamb))]

        elif cls.__approach == "MassConservationApproach":

            return [cls.__cons_laws_sympy_lamb[ii](*tuple([cls.__concentration_funs[j](*tuple(params_for_global_min)) for j in range(cls.__N)]))
                                 for ii in range(len(cls.__cons_laws_sympy_lamb))]

    @classmethod
    def run_direct_simulation(cls, params_for_global_min, parallel_flag, dir_path, itg, change_in_relative_error,
                              spec_index, left_multiplier, right_multiplier, comm, my_rank, num_cores, dir_sim_print_flag,
                              R, N, cons_laws_sympy_lamb, cons_laws_sympy, sympy_species, signal_index, signal, response,
                              ode_lambda_functions, jac_lambda_function, concentration_funs, approach):

        cls.__comm = comm
        cls.__my_rank = my_rank
        cls.__num_cores = num_cores
        cls.__dir_sim_print_flag = dir_sim_print_flag
        cls.__R = R
        cls.__N = N
        cls.__cons_laws_sympy_lamb = cons_laws_sympy_lamb
        cls.__signal_index = signal_index
        cls.__cons_laws_sympy = cons_laws_sympy
        cls.__sympy_species = sympy_species
        cls.__ode_lambda_functions = ode_lambda_functions
        cls.__jac_lambda_function = jac_lambda_function
        cls.__response = response
        cls.__signal = signal
        cls.__concentration_funs = concentration_funs
        cls.__approach = approach

        #TODO: import global modules needed, remember to do mpi_mod
        global itertools
        global sys
        global os
        global numpy
        global warnings
        global pandas
        global ggplot, aes, geom_line, ylim, scale_color_distiller, facet_wrap, theme_bw, geom_path, geom_point, labs, annotate

        import itertools
        import sys
        import os
        import numpy
        import warnings
        import pandas
        from plotnine import ggplot, aes, geom_line, ylim, scale_color_distiller, facet_wrap, theme_bw, geom_path, \
            geom_point, labs, annotate

        if parallel_flag is False and cls.__comm is None:

            # making the directory if it doesn't exist
            if not os.path.isdir(dir_path):
                os.mkdir(dir_path)

            print("Starting direct simulation ...")
            start_t = time.time()

        elif parallel_flag is True and cls.__comm is None:

            from mpi4py import MPI
            global mpi_mod
            from .mpi_routines import MPIRoutines as mpi_mod

            cls.__comm = MPI.COMM_WORLD
            cls.__my_rank = cls.__comm.Get_rank()
            cls.__num_cores = cls.__comm.Get_size()
            cls.__comm.Barrier()

            if not os.path.isdir(dir_path) and cls.__my_rank == 0:
                os.mkdir(dir_path)

            cls.__comm.Barrier()

            start_time = MPI.Wtime()

            if cls.__my_rank == 0:
                print("Starting direct simulation ...")
        elif cls.__comm is not None:

            from mpi4py import MPI
            global mpi_mod
            from .mpi_routines import MPIRoutines as mpi_mod

            if not os.path.isdir(dir_path) and cls.__my_rank == 0:
                os.mkdir(dir_path)

            cls.__comm.Barrier()

            params_for_global_min = cls.__comm.bcast(params_for_global_min, root=0)

            cls.__comm.Barrier()

            start_time = MPI.Wtime()

            if cls.__my_rank == 0:
                print("Starting direct simulation ...")

        else:
            print("Starting direct simulation ...")
            start_t = time.time()

        if len(params_for_global_min) == 0:
            print("The parameter sets provided has a length of zero, direct simulation cannot be ran.")
            sys.exit()

        viable_indices, viable_out_values, conservation_vals = cls.__find_viable_indices(params_for_global_min[0], itg,
                                                                                          spec_index,
                                                                                          change_in_relative_error)

        spec_index, fwd_scan_vals, rvrs_scan_vals, fwd_scan_index, rvrs_scan_index = cls.__initialize_direct_simulation(viable_indices, viable_out_values,
                                                                                                                        params_for_global_min[0], conservation_vals, itg,
                                                                                                                        change_in_relative_error, spec_index, left_multiplier,
                                                                                                                        right_multiplier)

        if cls.__dir_sim_print_flag:
            if cls.__comm is None:
                cls.__print_initial_conditions(fwd_scan_vals, rvrs_scan_vals)
            else:
                if cls.__my_rank == 0:
                    cls.__print_initial_conditions(fwd_scan_vals, rvrs_scan_vals)

        plot_flag = True

        for i in range(len(params_for_global_min)):

            if cls.__dir_sim_print_flag:
                if cls.__comm is None:
                    print(f"Conducting stability analysis of element {i} of the list provided ... ")
                else:
                    if cls.__my_rank == 0:
                        print(f"Conducting stability analysis of element {i} of the list provided ... ")

            conservation_vals = cls.construct_conservation_vals(params_for_global_min[i])
            # conservation_vals = [cls.__cons_laws_sympy_lamb[ii](*tuple(params_for_global_min[i][cls.__R:cls.__R + cls.__N]))
            #     for ii in range(len(cls.__cons_laws_sympy_lamb))]  # TODO: make general function

            con_law_value = conservation_vals[cls.__signal_index]
            change_left = con_law_value * left_multiplier
            change_right = con_law_value * right_multiplier
            pcp_scan = numpy.linspace(con_law_value - change_left, con_law_value + change_right, 100)

            forward_scan, reverse_scan = cls.__conduct_fwd_rvrs_scan(params_for_global_min[i], fwd_scan_vals,
                                                                     rvrs_scan_vals, pcp_scan, fwd_scan_index,
                                                                     rvrs_scan_index, spec_index, itg,
                                                                     change_in_relative_error)

            min_val, max_val = cls.__get_min_max_vals(pcp_scan, forward_scan, reverse_scan)

            count = 0
            scan_vals = pcp_scan
            while count < 5:
                if len([ii for ii in scan_vals if ii >= min_val and ii <= max_val]) < 10:

                    second_scan = numpy.linspace(min_val, max_val, 60)

                    forward_scan, reverse_scan = cls.__conduct_fwd_rvrs_scan(params_for_global_min[i], fwd_scan_vals,
                                                                              rvrs_scan_vals, second_scan, fwd_scan_index,
                                                                              rvrs_scan_index, spec_index, itg,
                                                                              change_in_relative_error)

                    min_val, max_val = cls.__get_min_max_vals(second_scan, forward_scan, reverse_scan)
                    scan_vals = second_scan
                    count += 1
                else:
                    break

            if count == 0:
                second_scan = pcp_scan

            if plot_flag:

                if cls.__comm is None:
                    cls.__plot_direct_simulation(second_scan, forward_scan, reverse_scan, dir_path, i)
                else:
                    if cls.__my_rank == 0:
                        cls.__plot_direct_simulation(second_scan, forward_scan, reverse_scan, dir_path, i)

        if cls.__comm is None:
            end_t = time.time()
            elapsed = end_t - start_t
            print("Elapsed time for direct simulation in seconds: " + str(elapsed))
        else:
            cls.__comm.Barrier()
            if cls.__my_rank == 0:
                end_time = MPI.Wtime()
                elapsed = end_time - start_time
                print(f"Elapsed time for direct simulation in seconds: {elapsed}")

    @classmethod
    def __find_viable_indices(cls, result_x, itg, spec_index, change_in_relative_error):

        conservation_vals = cls.construct_conservation_vals(result_x)
        # conservation_vals = [cls.__cons_laws_sympy_lamb[i](*tuple(result_x[cls.__R:cls.__R + cls.__N]))
        #                      for i in range(len(cls.__cons_laws_sympy_lamb))]                                    # TODO: make general function

        cons_laws_spec_info = []
        for i in cls.__cons_laws_sympy:
            spec_in_law = [ii for ii in cls.__sympy_species if ii in i.atoms()]
            spec_indices = [cls.__sympy_species.index(ii) for ii in spec_in_law]
            coeff_of_spec = [i.coeff(ii, 1) for ii in spec_in_law]
            cons_laws_spec_info.append([coeff_of_spec, spec_indices])

        temp_comb = list(itertools.product(*[i[1] for i in cons_laws_spec_info]))

        all_unique_comb = [temp_comb[i] for i in range(len(temp_comb)) if len(set(temp_comb[i])) == len(temp_comb[i])]

        viable_indices = []
        viable_out_values = []
        for i in all_unique_comb:
            initial_species_values = [0.0 for i in range(cls.__N)]

            for j in range(len(conservation_vals)):
                initial_species_values[i[j]] = conservation_vals[j]

            out = cls.__steady_state_finder(initial_species_values, result_x, spec_index, itg, change_in_relative_error)
            steady_cons = [cls.__cons_laws_sympy_lamb[i](*tuple(out)) for i in range(len(cls.__cons_laws_sympy_lamb))]

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

    @classmethod
    def __initialize_direct_simulation(cls, viable_indices, viable_out_values, result_x, conservation_vals, itg,
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

            fwd_scan_vals, rvrs_scan_vals, fwd_scan_index, rvrs_scan_index, fwd_ind, rvrs_ind = cls.__get_important_scan_vals(chosen_combo,
                                                                                                                              spec_index,
                                                                                                                              viable_out_values,
                                                                                                                              viable_indices)

            # determining if the forward or reverse scan is constant, if so, remove it as a viable combination
            con_law_value = conservation_vals[cls.__signal_index]
            # change = con_law_value*0.25
            # pcp_scan = numpy.linspace(con_law_value - change, con_law_value + change, 10)
            change_left = con_law_value * left_multiplier
            change_right = con_law_value * right_multiplier
            pcp_scan = numpy.linspace(con_law_value - change_left, con_law_value + change_right, 10)

            forward_scan, reverse_scan = cls.__conduct_fwd_rvrs_scan(result_x, fwd_scan_vals,
                                                                     rvrs_scan_vals, pcp_scan, fwd_scan_index,
                                                                     rvrs_scan_index, spec_index, itg,
                                                                     change_in_relative_error)

            combos, diff, stop_flag, combos_flag = cls.__get_new_combo(forward_scan, reverse_scan, fwd_ind, rvrs_ind, combos, diff)

            # if all combinations are thought to produce constant in time return the initial combo and continue
            if combos_flag:
                fwd_scan_vals, rvrs_scan_vals, fwd_scan_index, rvrs_scan_index, fwd_ind, rvrs_ind = cls.__get_important_scan_vals(
                    chosen_initial_combo,
                    spec_index,
                    viable_out_values,
                    viable_indices)
                break

        return spec_index, fwd_scan_vals, rvrs_scan_vals, fwd_scan_index, rvrs_scan_index

    @classmethod
    def __get_important_scan_vals(cls, chosen_combo, spec_index, viable_out_values, viable_indices):

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
        fwd_scan_index = [i[0] for i in fwd_scan_vals if i[1] == cls.__signal_index][0]
        # index to change in reverse scan
        rvrs_scan_index = [i[0] for i in rvrs_scan_vals if i[1] == cls.__signal_index][0]

        return fwd_scan_vals, rvrs_scan_vals, fwd_scan_index, rvrs_scan_index, fwd_ind, rvrs_ind

    @classmethod
    def __get_new_combo(cls, forward_scan, reverse_scan, fwd_ind, rvrs_ind, combos, diff):

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

    @classmethod
    def __get_min_max_vals(cls, pcp_scan, forward_scan, reverse_scan):

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

    @classmethod
    def __conduct_fwd_rvrs_scan(cls, result_x, fwd_scan_vals, rvrs_scan_vals, pcp_scan, fwd_scan_index,
                                rvrs_scan_index, spec_index, itg, change_in_relative_error):

        if cls.__comm is not None:
            pcp_scan = mpi_mod.distribute_list_of_points(pcp_scan, cls.__my_rank, cls.__num_cores, cls.__comm)

        conservation_vals = cls.construct_conservation_vals(result_x)
        # conservation_vals = [cls.__cons_laws_sympy_lamb[i](*tuple(result_x[cls.__R:cls.__R + cls.__N]))
        #                      for i in range(len(cls.__cons_laws_sympy_lamb))]                                    # TODO: make general function

        initial_species_values = [0.0 for i in range(cls.__N)]

        for i in fwd_scan_vals:
            initial_species_values[i[0]] = conservation_vals[i[1]]

        forward_scan = []
        for i in pcp_scan:
            initial_species_values[fwd_scan_index] = i
            steady_state = cls.__steady_state_finder(initial_species_values, result_x, spec_index, itg, change_in_relative_error)
            forward_scan.append(steady_state[spec_index])

        initial_species_values = [0.0 for i in range(cls.__N)]
        for i in rvrs_scan_vals:
            initial_species_values[i[0]] = conservation_vals[i[1]]

        reverse_scan = []
        for i in pcp_scan:
            initial_species_values[rvrs_scan_index] = i
            steady_state = cls.__steady_state_finder(initial_species_values, result_x, spec_index, itg, change_in_relative_error)
            reverse_scan.append(steady_state[spec_index])

        if cls.__comm is not None:
            list_forward_scan = mpi_mod.gather_list_of_values(forward_scan, cls.__comm, cls.__my_rank)
            list_reverse_scan = mpi_mod.gather_list_of_values(reverse_scan, cls.__comm, cls.__my_rank)

            list_forward_scan = cls.__comm.bcast(list_forward_scan, root=0)
            list_reverse_scan = cls.__comm.bcast(list_reverse_scan, root=0)

            cls.__comm.Barrier()

            return list_forward_scan, list_reverse_scan

        else:
            return forward_scan, reverse_scan

    @classmethod
    def __steady_state_finder(cls, initial_species_values, result_x, spec_index, itg, change_in_relative_error):

        def ff(t, cs, ks, ode_lambda_functions, jacobian):
            return [i(*tuple(ks), *tuple(cs)) for i in ode_lambda_functions]

        def jac_f(t, cs, ks, ode_lambda_functions, jacobian):
            return jacobian(*tuple(ks), *tuple(cs))

        len_time_interval = 100.0

        with numpy.errstate(divide='ignore', invalid='ignore'):
            out = itg.solve_ivp(ff, [0.0, len_time_interval], initial_species_values, args=(result_x[0:cls.__R], cls.__ode_lambda_functions, cls.__jac_lambda_function), jac=jac_f, method='BDF', rtol=1e-6, atol=1e-9, vectorized=True) #'RK45')  #'LSODA')
            y0 = out.y[:, -1]

        flag = True

        i = 1
        while flag:
            tspan = [0.0 + i*len_time_interval, len_time_interval + i*len_time_interval]
            try:
                with numpy.errstate(divide='ignore', invalid='ignore'):
                    out = itg.solve_ivp(ff, tspan, y0, args=(result_x[0:cls.__R], cls.__ode_lambda_functions, cls.__jac_lambda_function), jac=jac_f, method='BDF', rtol=1e-6, atol=1e-9, vectorized=True) #'RK45')  #'LSODA')
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

    @classmethod
    def __plot_direct_simulation(cls, pcp_scan, forward_scan, reverse_scan, path, index):

        out = pandas.DataFrame(columns=['dir', 'signal'] + [cls.__response])
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
             + aes(x='signal', y=cls.__response, color='dir')
             + labs(x=f"{cls.__signal} total", y=f"[{cls.__response}]", color="")
             + geom_path(size=2, alpha=0.5)
             + geom_point(color="black")
             + theme_bw()
             + geom_point(color="black"))
        g.save(filename=path + f"/sim_bif_diag_{index}.png", format="png", width=6, height=4, units='in', verbose=False)