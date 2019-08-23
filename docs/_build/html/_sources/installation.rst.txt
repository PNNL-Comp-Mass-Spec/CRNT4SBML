.. highlight:: shell

.. _my-installation-label:

============
Installation
============

Creating a Virtual Environment
--------------------------------

The preferred way to use CRNT4SBML is through a virtual environment. A virtual environment for Python is a self-contained
directory tree. This environment can have a particular version of Python and Python packages. This is very helpful as
it allows one to use different versions of Python and Python packages without their install conflicting with already
installed versions. Here we will give a brief description of creating a virtual environment for CRNT4SBML using
`virtualenv <https://virtualenv.pypa.io/en/latest/>`_. To begin we first obtain virtualenv through a `pip`_ install:

.. code-block:: console

    $ pip install virtualenv

Once virtualenv is installed, download the latest version of `Python 3.6 <https://www.python.org/downloads/>`_ (be sure
to take note of the download location). Next we will create a directory to hold all of the virtual environments that we
may create called "python\_environments":

.. code-block:: console

    $ mkdir python_environments

Now that we have virtualenv and Python version 3.6, we can create the virtual environment crnt4sbml\_env in the
directory python\_environments as follows:

.. code-block:: console

    $ cd python_environments
    $ virtualenv -p /path/to/python/3.6/interpreter crnt4sbml_env

The flag "-p" tells virtualenv to create an environment using a specific Python interpreter, note that if a standard
download of Python was followed, then the path to the interpreter is as follows "/usr/local/bin/python3.6". One can now
see a directory called "crnt4sbml\_env" is created in the directory python\_environments. We can now activate this
environment as follows:

.. code-block:: console

    $ source /path/to/python_environments/crnt4sbml_env/bin/activate

On the command line one should now see "(crnt4sbml_env)" on the left side of the command line, which indicates that one
is now working in the virtual environment. Once the environment is activated, one can now install CRNT4SBML as follows:

.. code-block:: console

    $ pip install crnt4sbml

note that this will install crnt4sbml in the virtual environment crnt4sbml_env. One can only use crnt4sbml within this
environment. If one wants to stop using the virtual environment, the following command can be used:

.. code-block:: console

    $ deactivate

"(base)" should show up on the left of the command line. One can then use the environment by using the "source" command
above.

Stable release
--------------

crnt4sbml can be obtained through a standard `pip`_ install as follows:

.. code-block:: console

    $ pip install crnt4sbml

This is the preferred method to install crnt4sbml, as it will always install the most recent stable release. Note that
crnt4sbml requires Python version 3.6.

.. _pip: https://pip.pypa.io

From sources
------------

The sources for crnt4sbml can be downloaded from the `Github repo`_.

You can also clone the public repository:

.. code-block:: console

    $ git clone https://github.com/PNNL-Comp-Mass-Spec/CRNT4SBML.git

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repo: https://github.com/PNNL-Comp-Mass-Spec/CRNT4SBML
