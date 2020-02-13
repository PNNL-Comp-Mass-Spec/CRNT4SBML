.. highlight:: shell

.. _parallel-crnt4sbml-label:

====================
Parallel CRNT4SBML
====================

Currently underdevelopment!!!!!!!!!!

Due to the nature of the optimization problem formed and the numerical continuation performed, some models can take
a long time to complete. In order to improve the user experience, we have developed parallel versions of the optimization
and numerical continuation routines for all approaches using `mpi4py <https://mpi4py.readthedocs.io/en/stable/>`_.

Installing the proper packages
+++++++++++++++++++++++++++++++++

Requirements for Parallel Version
-----------------------------------

- Python 3.7 (64-bit)
- networkx==2.3
- python-libsbml==5.18.0
- numpy==1.16.4
- sympy==1.4
- scipy==1.3.0
- matplotlib==3.1.0
- antimony==2.11.0
- rrplugins==1.2.2
- libroadrunner==1.5.2.1
- mpi4py==3.0.3

Creating a Virtual Environment
--------------------------------

The preferred way to use the parallel version of CRNT4SBML is through a virtual environment. First follow the steps
outlined in :ref:`my-installation-label` to create a virtual environment with the name mpi_crnt4sbml. Once this is done,
we can now activate this environment as follows:

On Mac:

.. code-block:: console

    $ source /path/to/python_environments/mpi_crnt4sbml/bin/activate

On Windows:

.. code-block:: console

    $ path\to\mpi_crnt4sbml\Scripts\activate

Note, in case you are using PowerShell, make sure its policy is updated by executing command as administrator
``Set-ExecutionPolicy RemoteSigned``. On the command line one should now see "(mpi_crnt4sbml)" on the left side of the
command line, which indicates that one is now working in the virtual environment.

Once the environment is activated, one can now install CRNT4SBML as follows:

.. code-block:: console

    $ pip install crnt4sbml

note that this will install crnt4sbml in the virtual environment mpi_crnt4sbml. One can only use crnt4sbml within this
environment.

One now needs to install mpi4py. Given mpi4py uses mpicc under the covers, we first need to install an MPI compiler
onto our system. This is done differently on MacOS and Windows.

On Mac:

    The simplest way to install mpicc on mac is to use homebrew. To begin, first install
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

On Windows:

    The simplest way to install a proper MPI compiler on Windows is to use Microsoft MPI. If not already installed,
    one should download Microsoft MPI version 10 or newer. At the time of creating this documentation, this could be
    done using the following `link <https://www.microsoft.com/en-us/download/details.aspx?id=57467>`_. Using the link
    click download and download msmpisetup.exe and run it. After the download, one should have a proper MPI compiler
    that is compatible with mpi4py.

    Note that for some users, one will also need to set the MSMPI path under User Variables. By default the Variable
    should be set to MSMPI_BIN and the Value should be ``C:\Program Files\Microsoft MPI\Bin``. This can be done
    following the instructions `here <https://www.computerhope.com/issues/ch000549.htm>`_.

    We are now ready to install mpi4py. With the virtual environment mpi_crnt4sbml activated, mpi4py can be installed
    as follows:

    .. code-block:: console

        $ pip install mpi4py


Parallel Mass Conservation Approach
+++++++++++++++++++++++++++++++++++++

To run the mass conservation approach create the following python script:

.. code-block:: python

   import crnt4sbml

   network = crnt4sbml.CRNT("path/to/Fig1Ci.xml")

   opt = network.get_mass_conservation_approach()

   bounds, concentration_bounds = opt.get_optimization_bounds()

   params_for_global_min, obj_fun_val_for_params, my_rank = opt.run_mpi_optimization(bounds=bounds,
                                                                                     concentration_bounds=concentration_bounds)

   multistable_param_ind, sample_points, plot_specifications = opt.run_mpi_greedy_continuity_analysis(species="s15", parameters=params_for_global_min,
                                                                                                      auto_parameters={'PrincipalContinuationParameter': 'C3'})

   opt.generate_report()

This will provide the following output along with creating the directory "num\_cont\_graphs" in your current
directory that contains multistability plots. Please note that runtimes and the number of multistability plots produced
may vary among different operating systems.

::

    Creating Equilibrium Manifold ...
    Creating Equilibrium Manifold ...
    Elapsed time for creating Equilibrium Manifold: 1.8794879999999994
    Elapsed time for creating Equilibrium Manifold: 1.8736319999999997

    Running feasible point method for 10 iterations ...
    Elapsed time for feasible point method: 0.788164

    Running the multistart optimization ...
    Elapsed time for multistart method in seconds: 3.1570869999999998

    Running continuity analysis ...
    Elapsed time for continuity analysis in seconds: 15.839016


    The number of feasible points that satisfy the constraints by core 1: 5
    Total feasible points that give F(x) = 0 by core 1: 2
    Total number of points that passed final_check by core 1: 2

    Number of multistability plots found by core 1: 0
    Elements in params_for_global_min that produce multistability found by core 1:
    []

    The number of feasible points that satisfy the constraints by core 0: 5
    Smallest value achieved by objective function: 0.0
    Total feasible points that give F(x) = 0 by core 0: 2
    Total number of points that passed final_check by core 0: 2

    Number of multistability plots found by core 0: 2
    Elements in params_for_global_min that produce multistability found by core 0:
    [0, 1]



Parallel Semi-diffusive Approach
+++++++++++++++++++++++++++++++++++++


Parallel General Approach
+++++++++++++++++++++++++++


.. _pip: https://pip.pypa.io

