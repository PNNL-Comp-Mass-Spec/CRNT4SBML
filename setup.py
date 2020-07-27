#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

setup_requirements = []
test_requirements = []

setup(
    author="Brandon Reyes",
    author_email='reyesb123@gmail.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.7',
    ],
    description="CRNT4SBML is an easily installable Python based package available on MacOS and Windows. CRNT4SBML is concentrated on providing a simple workflow for the testing of core CRNT methods directed at detecting bistability in cell signaling pathways endowed with mass action kinetics.",
    entry_points={
        'console_scripts': [
            'crnt4sbml=crnt4sbml.cli:main',
        ],
    },
    install_requires=['networkx==2.3',
                      'python-libsbml==5.18.0',
                      'numpy==1.16.4',
                      'sympy==1.4',
                      'scipy==1.4.1',
                      'matplotlib==3.1.1',
                      'plotnine==0.6.0'],
    extras_require={"Windows": ['antimony==2.11.0', 'libroadrunner==1.5.2.1', 'rrplugins==1.2.2'],
                    "MacOS": ['antimony==2.11.0', 'libroadrunner==1.5.2.1', 'rrplugins==1.2.2'],
                    "Linux": [''],
                    "MPIWindows": ['antimony==2.11.0', 'libroadrunner==1.5.2.1', 'rrplugins==1.2.2', 'mpi4py==3.0.3'],
                    "MPIMacOS": ['antimony==2.11.0', 'libroadrunner==1.5.2.1', 'rrplugins==1.2.2', 'mpi4py==3.0.3'],
                    "MPILinux": ['mpi4py==3.0.3']
                    },
    license="Apache Software License 2.0",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    platforms='ALL',
    keywords='crnt4sbml',
    name='crnt4sbml',
    packages=['crnt4sbml'],
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/PNNL-Comp-Mass-Spec/CRNT4SBML',
    version='0.0.12',
    zip_safe=False,
)
