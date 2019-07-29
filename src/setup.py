# Copyright (C) 2013-2015 Biao Li (biaol@bcm.edu) 
# License: GNU General Public License (http://www.gnu.org/licenses/)

from distutils.core import setup, Extension
try:
   from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
   from distutils.command.build_py import build_py

import sys, os
try:
    import argparse
except ImportError:
    sys.exit('This program requires Python 2.7.2 or higher, or Python 3.2 or higher. Please upgrade your version (%s) of Python and try again.' % (sys.version.split()[0]))

from src_simrareped import VERSION

setup(name = "RarePedSim",
    version = VERSION,
    description = "Simulation framework for generating family-based rare variant data",
    author = 'Biao Li and Gao Wang',
    author_email = 'biaol@bcm.edu',
    maintainer = 'Biao Li',
    maintainer_email = 'biaol@bcm.edu',
    py_modules = [
        'src_simrareped.__init__',
        'src_simrareped.simRarePed',
        'src_simrareped.utilityFunc',
        'src_simrareped.simulator',
        'src_simrareped.srv_batch_simrareped',
        'src_simrareped.gdata',
        'src_simrareped.test_spower',
	'src_simrareped.parallel',
    ],
    scripts = ['rarepedsim'],
    cmdclass = {'build_py': build_py },
    #package_dir = {'src_simrareped': 'simRarePed'},
    packages = ['src_simrareped'],
    package_data = {}
)
