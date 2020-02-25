"""
General basic routines for python interface to Fortran-Library
"""

''' Library '''
from typing import Any, List, Tuple
from pathlib import Path
from ctypes import *
import numpy # This also tells python to use numpy MKL
import numpy.random

''' Global variable '''
p_int    = POINTER(c_int   )
p_double = POINTER(c_double)

''' Auxiliary routine '''
def array2p(array:numpy.ndarray) -> Any:
    """ Return a C pointer converted from array """
    n = numpy.prod(array.shape)
    a = numpy.reshape(array, n)
    if a.dtype == int:
        p = cast((c_int*n)(), p_int)
    elif a.dtype == float:
        p = cast((c_double*n)(), p_double)
    for i in range(n): p[i] = a[i]
    return p

def p2array(p:Any, array:numpy.ndarray) -> None:
    """ Fill array with the value contained in p """
    n = numpy.prod(array.shape)
    a = numpy.reshape(array, n)
    for i in range(n): a[i] = p[i]
    array = numpy.reshape(a, array.shape)

''' Load libFL.so '''
FL = CDLL('libFL.so')