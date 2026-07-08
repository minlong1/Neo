"""
Shared import bootstrap for xps.py: pulls in input_arg's parsed CLI state
(including `timeing_mode`) via wildcard import, then times the remaining
third-party imports if timing mode is on. Import-time side effects are
intentional here (mirrors input_arg.py's own module-level argument
parsing) rather than something to refactor away mid-port.
"""

# Import Library
from PhysicsModules.XPS.xps_neo.input_arg import *
from PhysicsModules.XPS.xps_neo.helper import timecall

if timeing_mode:
# %matplotlib inline
    t1 = timecall()

from psutil import cpu_count
# Set the number of threads
import os
import sys
from operator import itemgetter
import numpy as np
import pathlib
import copy
import logging
import operator
import datetime
import random
import csv


import matplotlib as mpl
import matplotlib.pyplot as plt


if timeing_mode:
    initial_elapsed = timecall()- t1
    print('Inital import function took %.2f second' %initial_elapsed)
