# NEO

NEO is being refactored into a modular scientific computing framework. The new repository structure is organized around separate physics modules, 
with each module intended to contain its own source code, tests, example files, documentation, and setup instructions.

## Current Status

At the moment, the only supported module in this repository is EXAFS:

PhysicsModules/EXAFS

The EXAFS module contains the current working version of EXAFS Neo, which uses a genetic algorithm to fit Extended X-ray Absorption Fine Structure data.

Other modules may be added to the PhysicsModules directory in the future, but they are not currently supported unless explicitly documented.

## Repository Structure

The repository is organized around the following top-level layout:

NEO/
├── PhysicsModules/
│   └── EXAFS/
│       ├── tests/
│       ├── path_files/
│       ├── result/
│       └── README.md
├── setup.py
├── requirements.txt
└── README.md

The top-level README describes the overall NEO framework structure. Module-specific setup instructions, dependencies, example commands, and test information are maintained inside each module directory.

## Setting Up EXAFS

To install and run the EXAFS module, follow the README located inside the EXAFS module directory:

PhysicsModules/EXAFS/README.md

That README contains the currently supported setup process for EXAFS, including environment creation, dependency installation, and example usage.

A typical EXAFS test command from the repository root is:

exafs -i PhysicsModules/EXAFS/tests/cu_test_files/test_cu.ini