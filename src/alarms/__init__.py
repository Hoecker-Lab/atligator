"""This package contains reusable components of ALARMS within ATLIGATOR.

:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:date: 2022-01-19
"""

import warnings

warnings.warn("Alarms is depreciated. Please import the regarding modules from atligator directly.")
from atligator import chain_processing, construct_selection

# Just for keeping the imports
cp = chain_processing
cs = construct_selection
