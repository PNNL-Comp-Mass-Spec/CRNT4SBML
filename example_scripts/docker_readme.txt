Once you are in the directory where Dockerfile exists, do the following on 
the command line:
1. docker build -t crnt4sbml_image .
2. docker create --name crnt4sbml_container -t -i crnt4sbml_image /bin/bash

If you want to mount folders do the following for 2. :

2. docker create --name crnt4sbml_container --mount type=bind,source=/path/to/sbml_files,target=/home/crnt4sbml-user/sbml_files --mount type=bind,source=/host/path/to/example_scripts,target=/home/crnt4sbml-user/example_scripts -t -i crnt4sbml_image /bin/bash

To launch the container do the following: 

docker start -i crnt4sbml_container

Now that we are in the container, we can run any of the Python scripts for CRNT4SBML. 

Useful commands: 

       Show all containers:

       	    docker ps -a

       Show all images:

       	    docker images -a 
