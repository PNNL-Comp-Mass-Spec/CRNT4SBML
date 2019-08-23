===========
Quick Start
===========

To begin using CRNT4SBML, start by following the process outlined in :ref:`my-installation-label`. Once you have
correctly installed CRNT4SBML follow the steps below to obtain a general idea of how one can perform the mass conservation
and semi-diffusive approach of :cite:`irene`. If you are interested in running the Deficiency Zero and One theorems please
consult :ref:`my-basic-label`.

.. _quickstart-deficiency-label:

++++++++++++++++++++++++++++++++++++
Mass Conservation Approach Example
++++++++++++++++++++++++++++++++++++ 

In order to run the mass conservation approach one needs to first create an SBML file of the reaction network. The
SBML file representing the reaction network for this example is given by :download:`Fig1Ci.xml <../sbml_files/Fig1Ci.xml>`.
It is highly encouraged that the user consult :ref:`my-celldesigner-label` when considering their own individual network
as the format of the SBML file must follow a certain construction to be easily used by CRNT4SBML.

To run the mass conservation approach create the following python script:

.. code-block:: python

   import crnt4sbml

   network = crnt4sbml.CRNT("/path/to/Fig1Ci.xml")

   opt = network.get_mass_conservation_approach()

   bounds = [(1e-2, 1e2)]*12
   concentration_bounds = [(1e-2, 1e2)]*4

   params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds,
                                                                        concentration_bounds=concentration_bounds,
                                                                        iterations=15)

   multistable_param_ind = opt.run_greedy_continuity_analysis(species="s15", parameters=params_for_global_min,
                                                              auto_parameters={'PrincipalContinuationParameter': 'C3'})

   opt.generate_report()

This will provide the following output along with creating the directory "num\_cont\_graphs" in your current
directory that contains multistability plots. Please note that runtimes and the number of multistability plots produced
may vary among different operating systems. Please see :ref:`my-deficiency-label` for a more detailed explanation of
running the mass conservation approach and the provided output.

::

    Creating Equilibrium Manifold ...
    Elapsed time for creating Equilibrium Manifold: 0.294473

    Solving for species' concentrations ...
    Elapsed time for finding species' concentrations: 0.5590229999999998

    Running feasible point method for 15 iterations ...
    Elapsed time for feasible point method: 0.10936699999999977

    Running the multistart optimization ...

    Smallest value achieved by objective function: 0.0

    Elapsed time for multistart method: 17.663542

    Running continuity analysis ...
    Elapsed time for continuity analysis: 18.639788999999997

    The number of feasible points that satisfy the constraints: 15
    Total feasible points that give F(x) = 0: 6
    Total number of points that passed final_check: 6
    Number of multistability plots found: 6
    Elements in params_for_global_min that produce multistability:
    [0, 1, 2, 3, 4, 5]

.. _`quickstart-injectivity-label`:

+++++++++++++++++++++++++++++++++++++
Semi-diffusive Approach Example
+++++++++++++++++++++++++++++++++++++

To run the semi-diffusive approach one needs to create the SBML file specific for semi-diffusive networks. The SBML file
representing the reaction network for this example is given by :download:`Fig1Cii.xml <../sbml_files/Fig1Cii.xml>`. It
is highly encouraged that the user consult :ref:`my-celldesigner-label` when considering their own individual network as
the format of the SBML file must follow a certain construction to be easily used by crnt4sbml.

To run the semi-diffusive approach create the following python script:

.. code-block:: python

    import crnt4sbml

    network = crnt4sbml.CRNT("path/to/Fig1Cii.xml")

    opt = network.get_semi_diffusive_approach()

    bounds = [(1e-3, 1e2)]*12

    params_for_global_min, obj_fun_val_for_params = opt.run_optimization(bounds=bounds)

    multistable_param_ind = opt.run_greedy_continuity_analysis(species="s7", parameters=params_for_global_min,
                                                               auto_parameters={'PrincipalContinuationParameter': 're17'})

    opt.generate_report()


This will provide the following output along with creating the directory "num\_cont\_graphs" in your current
directory that contains multistability plots. Please note that runtimes and the number of multistability plots produced
may vary among different operating systems. Please see :ref:`my-injectivity-label` for a more detailed explanation of
running the semi-diffusive approach and the provided output.

::

    Running feasible point method for 10 iterations ...
    Elapsed time for feasible point method: 0.6484119999999995

    Running the multistart optimization ...

    Smallest value achieved by objective function: 0.0

    Elapsed time for multistart method: 43.25512

    Running continuity analysis ...
    Elapsed time for continuity analysis: 34.786898

    The number of feasible points that satisfy the constraints: 10
    Total feasible points that give F(x) = 0: 9
    Total number of points that passed final_check: 9
    Number of multistability plots found: 9
    Elements in params_for_global_min that produce multistability:
    [0, 1, 2, 3, 4, 5, 6, 7, 8]


