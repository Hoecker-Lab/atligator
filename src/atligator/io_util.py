"""This module contains utility functions for command line input and output.

:Author: Felix Schwaegerl <felix.schwaegerl@uni-bayreuth.de>
:Author: Josef Kynast <josef.kynast@uni-bayreuth.de>
:date: 2018-03-06
"""

import sys


def err(message: str) -> None:
    """Prints out an error message to stderr in a red color.
    :param message: The error message to print.
    """
    sys.stderr.write('\033[01;31m' + message + '\033[0;0m\n')


def reduce_float(floating_number: float, digits: int = 2):
    """
    Converts a floating number into a string and returns it with specified number of digits.
    :param floating_number: floating number to use
    :param digits: digits after the decimal point to display
    :return: string of floating number with 'digits' digits after decimal point
    """
    return ("{0:." + str(digits) + "f}").format(floating_number)


class Colors:
    """
    Class for coloring printed text. When initialized as colors allows to
    - adjust color/style of printed text by 'print(color.RED, "hello")'
    - change all prints by color.red (equal to 'print(color.RED, end="")')
    - change back to normal by color.white or color.WHITE (as a string in print)
    """
    LMAGENTA: str = '\033[95m'
    BLUE: str = '\033[94m'
    GREEN: str = '\033[32m'
    LGREEN: str = '\033[92m'
    YELLOW: str = '\033[93m'
    LRED: str = '\033[91m'
    GREY: str = '\033[90m'
    RED: str = '\033[31m'
    BLACK: str = '\033[30m'
    WHITE: str = '\033[0m'
    BOLD: str = '\033[1m'
    UNDERLINE: str = '\033[4m'
    DIM: str = '\033[2m'
    BLINK: str = '\033[5m'
    HIDDEN: str = '\033[8m'
    BACKGROUND: str = '\033[7m'

    def pc(self, item) -> None:
        """
        Prints the attribute of choice.
        :param item: attribute to print (e.g. a color)
        :return: None
        """
        print(getattr(self, item.upper()), end="")


class Blank(object):
    """
    Class as a blank placeholder for importing Colors class without changing colors or style.
    Use similar to: 'from atligator.io_util import Blank as Colors'  When coloring option is disabled.
    """
    def __getattribute__(self, item):
        return ''


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = None

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self._original_stdout
