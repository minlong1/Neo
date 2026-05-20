# EXAFS Neo

#### Versions: 2.0.14

#### Last update: Feb 4, 2026

<!-- ![example workflow](https://github.com/laumiulun/EXAFS_Neo/actions/workflows/<WORKFLOW_FILE>/badge.svg) -->

[![Test with Ubuntu, Mamba](https://github.com/laumiulun/EXAFS_Neo/actions/workflows/test_ubuntu.yml/badge.svg?branch=devel)](https://github.com/laumiulun/EXAFS_Neo/actions/workflows/test_ubuntu.yml)[![Test with Windows, Mamba](https://github.com/laumiulun/EXAFS_Neo/actions/workflows/test_windows.yml/badge.svg?branch=devel)](https://github.com/laumiulun/EXAFS_Neo/actions/workflows/test_windows.yml)


EXAFS Neo utilize Genetic algorithm in fitting Extended X-ray absorption fine structure(EXAFS).

## Prerequisites

It is highly recommend to utilize `anaconda` or `pip` package managers to prevent unforeseen dependency conflicts. EXAFS Neo uses [`larch`](https://xraypy.github.io/xraylarch/) to process the x-ray spectrum.

- Python: => 3.9
- Numpy: => 1.20
- Scipy: => 1.6
- Larch: > 0.9.47
- Matplotlib: > 3.0

It is highly recommend to create a new environment in `anaconda` to run EXAFS Neo to prevent packages conflicts. For `Windows` operating system, if you encounter an issue requiring "Microsoft C++ 14.0 or greater is needed", please download the tools at the following location [`C++ Tools`](https://visualstudio.microsoft.com/visual-cpp-build-tools/) and make sure to select C++ build tools during installation process.

if you are on a Mac (either Intel or M1), you need to make sure that xcode command line tools is install, if not input this command into terminal:

        xcode-select --install

## Dependencies
<details>
<summary>Notes: This is usually not needed unless you are trying to install via source</summary>

EXAFS Neo requires the following dependencies to run:

        # Create new anaconda environment
        conda create -y --name exafs python=>3.9.10
        conda activate exafs
        conda install -y -c conda-forge "numpy>=1.23" "scipy>=1.8" "matplotlib>=3.6" "h5py>=3.5" "wxpython>=4.1" scikit-image scikit-learn pycifrw pandas jupyter plotly pyparsing pytest pytest-cov coverage
        pip install lmfit peakutils pyepics pyshortcuts termcolor sphinx dill pycifrw xraydb wxmplot wxutils fabio silx imageio charset-normalizer
        pip install xraylarch

</details>

## Installations

### Linux/Mac

**Conda** 

        conda create --name "exafs_neo" "python=3.12"
        conda activate exafs_neo
        pip install exafs_neo

**PyPI**

        python -m venv .venv
        source .venv/bin/activate
        pip install exafs_neo

<details>
<summary>

### Windows Installation (WSL)
</summary>

For Windows users, we recommend using the Windows Subsystem for Linux (WSL) to ensure full compatibility with GUI applications.

Note: Update Graphics Drivers: Follow the official guide to run Linux GUI apps on the Windows Subsystem for Linux.

1. Install WSL: Open PowerShell as Administrator and run:

    ```Powershell
    wsl --install
    ```
2. Install Visual Studio: Ensure you have Visual Studio installed (often required for compiling C++ extensions).
        
3. Install System Dependencies: Open your WSL terminal (Ubuntu/Debian) and run the following commands:

    ```Bash
    sudo apt update
    sudo apt install tk
    sudo apt install python3-tk
    sudo apt install pipx
    ```
4. Install EXAFS Neo:

    ```Bash
    sudo pipx install exafs-neo
    ```
</details>

### Source
To install EXAFS Neo, simply clone the repo:

        git clone https://github.com/laumiulun/EXAFS_Neo.git
        cd EXAFS_Neo/
        pip install .

## Usage

To run the included test suite, use the following command:

        ./run_tests

To run a sample test, make sure the environment is set correctly, and select an input file:

         exafs_neo -i test_inputs/test_temp.ini

Alternatively, you can also run EXAFS Neo in a jupyter notebook, please follow the example in the example/jupyter folder

## GUI

We also have provided a GUI for use in additions to our program, with additional helper script to facilitate post-analysis. To use the GUI:

        exafs_neo_gui

The datafile requires header contain at least either a combination of (k, chi) or (energy, mu). It also requires a minimum of one newline for it to work correctly. An example of the correct header is as follows:

        #---------------------------------------------------------------------
        #  k chi chik chik2 chik3 win energy

## Self adsorption correction

EXAFS also provides an internal option to perform self-adsorption on the sample file using Booth et al. correction. This is performed using git submodules:

        git submodule update --init --recursive
        cd contrib/sabcor/
        make

## Update

EXAFS Neo is under active development, to update the code after pulling from the repository:

        git pull --rebase
        python setup.py install



## Video Demonstrations

You can see a list of video demonstrations of the EXAFS Neo package presented, future presentation related to this software will be posted as they are available

<!-- - https://www.youtube.com/playlist?list=PLqZCvArs4yF8IrREQ3AzZJX2N-IRAPEmy [Aug 23, 2021] (IIT EXAFS Workshop 2021)
- https://youtu.be/KwhItvwhapg [Feb 15, 2021] (University of Washington)
- https://youtu.be/jqISqq_FFR8 [Dec 10, 2020] (Canadian Light Source) -->

- [IIT EXAFS Workshop 2021](https://www.youtube.com/playlist?list=PLqZCvArs4yF8IrREQ3AzZJX2N-IRAPEmy) (Aug 23, 2021)
- [University of Washington](https://youtu.be/KwhItvwhapg) (Feb 15, 2021)
- [Canadian Light Source](https://youtu.be/jqISqq_FFR8) (Dec 10, 2020)

## Citation

Jeff Terry, Miu Lun Lau, Jiateng Sun, Chang Xu, Bryan Hendricks, Julia Kise, Mrinalini Lnu, Sanchayni Bagade, Shail Shah, Priyanka Makhijani, Adithya Karantha, Travis Boltz, Max Oellien, Matthew Adas, Shlomo Argamon, Min Long, and Donna Post Guillen, “Analysis of Extended X-ray Absorption Fine Structure (EXAFS) Data Using Artificial Intelligence Techniques,” Applied Surface Science 547, 149059 <https://doi.org/10.1016/j.apsusc.2021.149059> (2021).

    @article{terry2021analysis,
      title={Analysis of extended X-ray absorption fine structure (EXAFS) data using artificial intelligence techniques},
      author={Terry, Jeff and Lau, Miu Lun and Sun, Jiateng and Xu, Chang and Hendricks, Bryan and Kise, Julia and Lnu, Mrinalini and Bagade, Sanchayni and Shah, Shail and Makhijani, Priyanka and others},
      journal={Applied Surface Science},
      volume={547},
      pages={149059},
      year={2021},
      publisher={Elsevier}
    }
