import os
import pickle
import numpy
import antimony
import roadrunner
import rrplugins
import sys

roadrunner.Logger.setLevel(roadrunner.Logger.LOG_ERROR)
roadrunner.Logger.disableLogging()
roadrunner.Logger.disableConsoleLogging()
roadrunner.Logger.disableFileLogging()
rrplugins.setLogLevel('error')

stderr_fileno = sys.stderr.fileno()
stderr_save = os.dup(stderr_fileno)
stderr_pipe = os.pipe()
os.dup2(stderr_pipe[1], stderr_fileno)
os.close(stderr_pipe[1])


# functions taken from Tellurium!! Give them
# credit, they deserve it!
#################################################
def __check_antimony_return_code(code):
    if code < 0:
        raise Exception('Antimony: {}'.format(antimony.getLastError()))


def __antimony_to_sbml(ant):
    try:
        isfile = os.path.isfile(ant)
    except ValueError:
        isfile = False
    if isfile:
        code = antimony.loadAntimonyFile(ant)
    else:
        code = antimony.loadAntimonyString(ant)
    __check_antimony_return_code(code)
    mid = antimony.getMainModuleName()
    return antimony.getSBMLString(mid)


def __loada(ant):
    return __load_antimony_model(ant)


def __load_antimony_model(ant):
    sbml = __antimony_to_sbml(ant)
    return roadrunner.RoadRunner(sbml)


with open('input_arguments.pickle', 'rb') as pickle_file:
    input_arguments = pickle.loads(pickle_file.read())

ant_str = input_arguments[0]
direction = input_arguments[1]
auto = rrplugins.Plugin("tel_auto2000")
auto_parameters = input_arguments[2]

antimony_r = __loada(ant_str)

# # making the directory auto_fort_files if is does not exist
# if not os.path.isdir("./auto_fort_files"):
#     os.mkdir("./auto_fort_files")

auto.setProperty("SBML", antimony_r.getCurrentSBML())
auto.setProperty("ScanDirection", direction)
auto.setProperty("PreSimulation", "True")
auto.setProperty("PreSimulationDuration", 1.0)
auto.setProperty('KeepTempFiles', True)
auto.setProperty("TempFolder", "auto_fort_files")

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

output_arguments = [pts, lbls, ant_float_ids, flag]

if os.path.exists("output_arguments.pickle"):
    os.remove("output_arguments.pickle")
    with open('output_arguments.pickle', 'wb') as outf:
        outf.write(pickle.dumps(output_arguments))
else:
    with open('output_arguments.pickle', 'wb') as outf:
        outf.write(pickle.dumps(output_arguments))

os.close(stderr_pipe[0])
os.dup2(stderr_save, stderr_fileno)
os.close(stderr_save)
os.close(stderr_fileno)
