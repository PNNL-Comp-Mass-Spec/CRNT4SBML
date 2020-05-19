FROM python:3.7.4-slim

RUN DEBIAN_FRONTEND=noninteractive apt-get update \
&& apt-get install -y libopenmpi-dev \
&& apt-get install -y git \
&& pip install git+https://github.com/PNNL-Comp-Mass-Spec/CRNT4SBML.git#egg=crnt4sbml[Linux] \
&& rm -rf /var/lib/apt \
&& useradd -ms /bin/bash crnt4sbml-user 

USER crnt4sbml-user 
WORKDIR /home/crnt4sbml-user

