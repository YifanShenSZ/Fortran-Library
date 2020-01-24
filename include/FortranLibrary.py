"""
Python interface to Fortran-Library
Note that python cannot change function argument value by passing reference,
so such fortran functions now return output values in python style
"""

''' Library '''
import numpy
from ctypes import *

''' Global variable '''
p_int    = POINTER(c_int   )
p_double = POINTER(c_double)

''' Load libFL.so '''
FL = CDLL('libFL.so')

''' Rewrap fortran functions '''
# There are 2 reasons:
#     1. Python recognizes fortran functions by their compiled names:
#        fortran function Func contained in module Mod will be renamed as mod_mp_func_
#        ( You may use command nm to view the compiled names: nm libFL.so )
#     2. Most python function arguments are immutable:
#        such argument value cannot be changed by passing reference
# So we rename fortran functions back to their origin,
# and return output values in python style

# General
ShowTime = FL.general_mp_showtime_

dScientificNotation_temp = FL.general_mp_dscientificnotation_
dScientificNotation_temp.argtypes = [p_double, p_int]
def dScientificNotation(x: float, i: int):
    x_copy = c_double(x); i_copy = c_int(i)
    px = p_double(x_copy); pi = p_int(i_copy)
    dScientificNotation_temp(px,pi)
    return px[0],pi[0]

# FortranLibrary
TestFortranLibrary = FL.fortranlibrary_mp_testfortranlibrary_
