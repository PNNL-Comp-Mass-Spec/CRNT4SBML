.. highlight:: shell

.. _parallel-crnt4sbml-label:

====================
Parallel CRNT4SBML
====================

Due to the nature of the optimization problem formed, some models can take a long time to complete. In order to improve
the user experience, we have developed parallel versions of the optimization routine for all approaches
using `mpi4py <https://mpi4py.readthedocs.io/en/stable/>`_.

Installing the proper packages
+++++++++++++++++++++++++++++++++

Base Requirements for Parallel Version
-----------------------------------------

- Python 3.7 (64-bit)
- networkx==2.3
- python-libsbml==5.18.0
- numpy==1.16.4
- sympy==1.4
- scipy==1.4.1
- matplotlib==3.1.1
- plotnine==0.6.0
- mpi4py==3.0.3

MacOS and Windows
++++++++++++++++++++
- antimony==2.11.0
- rrplugins==1.2.2
- libroadrunner==1.5.2.1

Creating a Virtual Environment
--------------------------------

The preferred way to use the parallel version of CRNT4SBML is through a virtual environment. First follow the steps
outlined in :ref:`my-installation-label` to create a virtual environment with the name mpi_crnt4sbml. Once this is done,
we can now activate this environment as follows:

On MacOS and Linux:

.. code-block:: console

    $ source /path/to/python_environments/mpi_crnt4sbml/bin/activate

On Windows:

.. code-block:: console

    $ path\to\mpi_crnt4sbml\Scripts\activate

Note, in case you are using PowerShell, make sure its policy is updated by executing command as administrator
``Set-ExecutionPolicy RemoteSigned``. On the command line one should now see "(mpi_crnt4sbml)" on the left side of the
command line, which indicates that one is now working in the virtual environment.

One now needs to install mpi4py. Given mpi4py uses mpicc under the covers, we first need to install an MPI compiler
onto our system. This is done differently on MacOS, Linux, and Windows.

On MacOS:

    The simplest way to install mpicc on MacOS is to use homebrew. To begin, first install
    `homebrew <https://docs.brew.sh/Installation>`_. Then, we need to install open-mpi. This is done in the terminal
    as follows:

    .. code-block:: console

        $ brew install open-mpi

    Be sure to take note of the install location of open-mpi. We now need to set the environment variable for the MPI
    compiler. This is done as follows in the terminal (take note that here we are using version 4.0.2 of open-mpi):

    .. code-block:: console

        $ export MPICC=path/to/open-mpi/4.0.2/bin/mpicc

    If a standard install was followed, "path/to/" can be replaced with "/usr/local/Cellar/". We are now ready to
    install mpi4py. With the virtual environment mpi_crnt4sbml activated, mpi4py can be installed as follows:

    .. code-block:: console

        $ pip install mpi4py

On Linux:

    The simplest way to install an MPI compiler on Linux is to install open-mpi. This is done in the terminal as
    follows (note that one may need to use sudo):

    .. code-block:: console

        $ apt-get install -y libopenmpi-dev

On Windows:

    The simplest way to install a proper MPI compiler on Windows is to use Microsoft MPI. If not already installed,
    one should download Microsoft MPI version 10 or newer. At the time of creating this documentation, this could be
    done using the following `link <https://www.microsoft.com/en-us/download/details.aspx?id=57467>`_. Using the link
    click download and download msmpisetup.exe and run it. After the download, one should have a proper MPI compiler
    that is compatible with mpi4py.

    Note that for some users, one will also need to set the MSMPI path under User Variables. By default the Variable
    should be set to MSMPI_BIN and the Value should be ``C:\Program Files\Microsoft MPI\Bin``. This can be done
    following the instructions `here <https://www.computerhope.com/issues/ch000549.htm>`_.

Once the environment is activated, one can now install a parallel CRNT4SBML as follows:

On MacOS:
    .. code-block:: console

        $ pip install crnt4sbml[MPIMacOS]

On Windows:
    .. code-block:: console

        $ pip install crnt4sbml[MPIWindows]

On Linux (numerical continuation is unavailable for Linux):
    .. code-block:: console

        $ pip install crnt4sbml[MPILinux]

note that this will install crnt4sbml in the virtual environment mpi_crnt4sbml. One can only use crnt4sbml within this
environment.

Parallel Mass Conservation Approach
+++++++++++++++++++++++++++++++++++++

To run the optimization for the mass conservation approach create the following python script named mpi_run.py:

.. code-block:: python

   import crnt4sbml
   import numpy

   network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")

   approach = network.get_mass_conservation_approach()

   bounds, concentration_bounds = approach.get_optimization_bounds()

   params_for_global_min, obj_fun_val_for_params = approach.run_optimization(bounds=bounds, concentration_bounds=concentration_bounds,
                                                                             parallel_flag=True)

   if approach.get_my_rank() == 0:
       numpy.save('params.npy', params_for_global_min)

   approach.generate_report()

You can then run the script from the console using 2 cores using the following command:

.. code-block:: console

    $ mpiexec -np 2 python mpi_run.py


This will provide the following output along with saving the params_for_global_min to the file params.npy in the current
directory. You can then load in params.npy and run a serial version of the numerical continuation. Please note that
runtimes may vary among different operating systems.

::

    Creating Equilibrium Manifold ...
    Creating Equilibrium Manifold ...
    Elapsed time for creating Equilibrium Manifold: 2.06032
    Elapsed time for creating Equilibrium Manifold: 2.0805279999999993

    Running feasible point method for 10 iterations ...
    Elapsed time for feasible point method: 1.024346

    Running the multistart optimization method ...
    Elapsed time for multistart method: 3.5696950000000003

    Smallest value achieved by objective function: 0.0
    4 point(s) passed the optimization criteria.

Parallel Semi-diffusive Approach
+++++++++++++++++++++++++++++++++++++

To run the optimization for the semi-diffusive approach create the following python script named mpi_run.py:

.. code-block:: python

   import crnt4sbml
   import numpy

   network = crnt4sbml.CRNT("path/to/Fig1Cii.xml")

   approach = network.get_semi_diffusive_approach()

   bounds = approach.get_optimization_bounds()

   params_for_global_min, obj_fun_val_for_params = approach.run_optimization(bounds=bounds, parallel_flag=True)

   if approach.get_my_rank() == 0:
       numpy.save('params.npy', params_for_global_min)

   approach.generate_report()

You can then run the script from the console using 2 cores using the following command:

.. code-block:: console

    $ mpiexec -np 2 python mpi_run.py

This will provide the following output along with saving the params_for_global_min to the file params.npy in the current
directory. You can then load in params.npy and run a serial version of the numerical continuation. Please note that
runtimes may vary among different operating systems.

::

    Running feasible point method for 10 iterations ...
    Elapsed time for feasible point method: 0.38841

    Running the multistart optimization method ...
    Elapsed time for multistart method: 17.330986000000003

    Smallest value achieved by objective function: 0.0
    9 point(s) passed the optimization criteria.

.. _parallel-gen-app-label:

Parallel General Approach
+++++++++++++++++++++++++++

Further libraries required
---------------------------

- plotnine==0.6.0

To run the optimization and direct simulation bistability anaylsis for the general approach create the following
python script named mpi_run.py:

.. code-block:: python

   import crnt4sbml

   network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")

   approach = network.get_general_approach()

   bnds = approach.get_optimization_bounds()

   approach.initialize_general_approach(signal="C3", response="s15", fix_reactions=True)

   params_for_global_min, obj_fun_vals = approach.run_optimization(bounds=bnds, dual_annealing_iters=100, confidence_level_flag=True,
                                                                   parallel_flag=True)

   approach.run_direct_simulation(params_for_global_min, parallel_flag=True)

   approach.generate_report()

You can then run the script from the console using 4 cores using the following command:

.. code-block:: console

    $ mpiexec -np 4 python mpi_run.py

This will provide the following output along with saving the direct simulation plots in the directory path
./dir_sim_graphs. Please note that runtimes may vary among different operating systems.

::

    Running the multistart optimization method ...
    Elapsed time for multistart method: 10.842817

    Starting direct simulation ...
    Elapsed time for direct simulation in seconds: 270.852905
    It was found that 0.0 is the minimum objective function value with a confidence level of 1.0 .
    9 point(s) passed the optimization criteria.

.. _pip: https://pip.pypa.io

