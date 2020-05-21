
.. highlight:: shell

.. _docker-docs-label:

======================
Docker and CRNT4SBML
======================

To further the accessibility of CRNT4SBML, we have created a Dockerfile for CRNT4SBML. This allows one to use the full
Linux version of CRNT4SBML. `Docker <https://www.docker.com/why-docker>`_ is a software platform that uses OS-level
virtualization to deliver software in packages called containers. Although there are many reasons to use Docker, our main
use case will be to provide our users with a simple install of CRNT4SBML. To begin, first install Docker and then download
:download:`Dockerfile <../example_scripts/Dockerfile>` into the directory of your choice.

Once you are in the directory where Dockerfile exists, one can create an image of CRNT4SBML named "crnt4sbml_image" by
completing the following in a terminal:

.. code-block:: console

    $ docker build -t crnt4sbml_image .

Using this image, we can then create a basic container named "crnt4sbml_container" using the following command:

.. code-block:: console

    $ docker create --name crnt4sbml_container -t -i crnt4sbml_image /bin/bash

Alternatively, if one would like to mount the folders "sbml_files" and "example_scripts" of the host machine to the
container upon creation one can do the following:

.. code-block:: console

    $ docker create --name crnt4sbml_container --mount type=bind,source=/path/to/sbml_files,target=/home/crnt4sbml-user/sbml_files --mount type=bind,source=/host/path/to/example_scripts,target=/home/crnt4sbml-user/example_scripts -t -i crnt4sbml_image /bin/bash

This will allow the user to easily access and add both sbml files and python scripts between Docker and the host machine.
To launch the container do the following:

.. code-block:: console

    $ docker start -i crnt4sbml_container

Now that we are in the container, we can run any of the Python scripts for CRNT4SBML that are available for Linux
(in particular :func:`crnt4sbml.GeneralApproach` and :func:`crnt4sbml.MassConservationApproach` without numerical
continuation).

Useful commands:

    Show all containers:

    .. code-block:: console

        $ docker ps -a

    Show all images:

    .. code-block:: console

        $ docker images -a