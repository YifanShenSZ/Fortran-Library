"""
Python interface to Fortran-Library
Note that some python data types are immutable,
so such fortran routines now return output values in python style
"""

''' Library '''
from typing import Any, List, Tuple
from ctypes import *
import numpy # This tells python to use MKL in numpy

''' Global variable '''
p_int    = POINTER(c_int   )
p_double = POINTER(c_double)

''' Auxiliary routine '''
# Return a C pointer converted from array
def array2p(array:numpy.ndarray) -> Any:
    n = numpy.prod(array.shape)
    a = numpy.reshape(array, n)
    if a.dtype == int:
        p = cast((c_int*n)(), p_int)
    elif a.dtype == float:
        p = cast((c_double*n)(), p_double)
    for i in range(n): p[i] = a[i]
    return p

# Fill array with the value contained in p
def p2array(p:Any, array:numpy.ndarray) -> None:
    n = numpy.prod(array.shape)
    a = numpy.reshape(array, n)
    for i in range(n): a[i] = p[i]
    array = numpy.reshape(a, array.shape)

''' Load libFL.so '''
FL = CDLL('libFL.so')

''' Rewrap fortran functions '''
# There are 2 reasons:
#     1. Python recognizes fortran functions by their compiled names:
#        fortran function Func contained in module Mod will be renamed as mod_mp_func_
#        ( You may use command nm to view the compiled names: nm libFL.so )
#     2. Some python data types are immutable:
#        such argument value cannot be changed by passing reference
# So we rename fortran functions back to their origin,
# and return immutable output values in python style

""" General """
ShowTime = FL.general_mp_showtime_

def dScientificNotation(x:float, i:int) -> (float, int):
    x_copy = c_double(x); i_copy = c_int(i)
    FL.general_mp_dscientificnotation_(byref(x_copy),byref(i_copy))
    return x_copy.value, i_copy.value

""" Mathematics """
# Fail to fetch constant:
# fortran parameter no longer has explicit name in .so
# So we directly provide them in python
AMUInAU  = 1822.888486192
AInAU    = 1.8897261339212517
cm_1InAu = 0.000004556335830019422

""" Chemistry """
def Avogadro_Vibration(NAtoms:int, symbol:numpy.ndarray, r:numpy.ndarray, vibdim:int, freq:numpy.ndarray, mode:numpy.ndarray, FileName='avogadro.log') -> None:
    p_symbol = cast(((c_char*2)*NAtoms)(), POINTER(c_char*2))
    for i in range(NAtoms):
        p_symbol[i][0] = (symbol[i][0]).encode('ascii')
        if len(symbol[i]) < 2:
            p_symbol[i][1] = ' '.encode('ascii')
        else:
            p_symbol[i][1] = symbol[i][1].encode('ascii')
    p_r = array2p(r); p_freq = array2p(freq); p_mode = array2p(mode)
    n = len(FileName)
    f = (c_char*n)(); f.value = FileName.encode('ascii')
    FL.chemistry_mp_avogadro_vibration_\
        (byref(c_int(NAtoms)), p_symbol, p_r, byref(c_int(vibdim)), p_freq, p_mode, f, n)

""" GeometryTransformation """
# Fail to fetch type InvolvedMotion:
# InvolvedMotion.atom is an intel fortran array descriptor without c coutner part
# So we provide a similar python class
class InvolvedMotion:
    def __init__(self, MotionType:str, coeff:float, atom:List):
        self.type  = MotionType
        self.coeff = coeff
        self.atom  = atom.copy()

# Fail to fetch type InternalCoordinateDefinition:
# It depends on type InvolvedMotion
# So we provide a similar python class
class InternalCoordinateDefinition:
    def __init__(self):
        # We no longer provide self.NMotions as it = len(self.motion)
        self.motion = []

# ---------- Cartesian <-> Internal ----------

# Fail to fetch GeometryTransformation_IntCDef:
# It depends on class InternalCoordinateDefinition
# So we provide a function to load internal coordinate definition in python
def FetchInternalCoordinateDefinition(definition:str) -> (int, List):
    intdim=0; intcdef = []
    if definition == 'Columbus7':
        # First line is always 'TEXAS'
        # New internal coordinate line starts with 'K'
        with open('intcfl','r') as f: lines=f.readlines()
        for i in range(1,len(lines)):
            if lines[i][0] == 'K':
                intdim += 1
                intcdef.append(InternalCoordinateDefinition())
            if lines[i][20:24] == 'STRE':
                MotionType = 'stretching'
                atom = [int(lines[i][28:33]), int(lines[i][34:43])]
            elif lines[i][20:24] == 'BEND':
                MotionType = 'bending'
                atom = [int(lines[i][28:34]), int(lines[i][45:54]), int(lines[i][35:44])]
            elif lines[i][20:24] == 'TORS':
                MotionType = 'torsion'
                atom = [int(lines[i][28:34]), int(lines[i][35:44]), int(lines[i][45:54]), int(lines[i][55:64])]
            elif lines[i][20:23] == 'OUT':
                MotionType = 'OutOfPlane'
                atom = [int(lines[i][28:34]), int(lines[i][45:54]), int(lines[i][55:64]), int(lines[i][35:44])]
            if lines[i][10:20] == '          ': coeff = 1.0
            else: coeff = float(lines[i][10:20])
            intcdef[intdim-1].motion.append(InvolvedMotion(MotionType,coeff,atom))
    else:
        # First 6 spaces of a line are reserved to indicate the start of new internal coordinate
        # Example:
        #  coor |   coeff   |    type     |      atom
        # --------------------------------------------------
        #      1    1.000000    stretching     1     2          # Comment
        #           1.000000    stretching     1     3
        #      2    1.000000    stretching     1     2
        #          -1.000000    stretching     1     3
        #      3    1.000000       bending     2     1     3
        with open('InternalCoordinateDefinition','r') as f: lines=f.readlines()
    return intdim, intcdef

def DefineInternalCoordinate(definition:str) -> int:
    n = len(definition)
    f = (c_char*n)(); f.value = definition.encode('ascii')
    intdim = FL.geometrytransformation_mp_defineinternalcoordinate_(f,n)
    return intdim

# ========== Cartesian -> Internal ==========

def InternalCoordinateq(r:numpy.ndarray, q:numpy.ndarray, cartdim:int, intdim:int) -> None:
    p_r = array2p(r)
    p_q = array2p(q)
    FL.geometrytransformation_mp_internalcoordinateq_\
        (p_r, p_q, byref(c_int(cartdim)), byref(c_int(intdim)))
    p2array(p_q, q)

# Due to row- and column-major difference, python fetchs B^T
def WilsonBMatrixAndInternalCoordinateq(r:numpy.ndarray, BT:numpy.ndarray, q:numpy.ndarray, cartdim:int, intdim:int) -> None:
    p_r = array2p(r)
    p_BT = array2p(BT); p_q = array2p(q)
    FL.geometrytransformation_mp_wilsonbmatrixandinternalcoordinateq_\
        (p_r, p_BT, p_q, byref(c_int(cartdim)), byref(c_int(intdim)))
    p2array(p_BT, BT); p2array(p_q, q)

# =================== End ===================

# ------------------- End --------------------

# --------------- Normal mode ----------------
# Due to row- and column-major difference, python
#     throws: H^T (H^T = H), B^T
#     fetchs: L^T, (L^-1)^T
def WilsonGFMethod(H:numpy.ndarray, BT:numpy.ndarray, mass:numpy.ndarray, freq:numpy.ndarray, LT:numpy.ndarray, LinvT:numpy.ndarray, intdim:int, NAtoms:int) -> None:
    p_H = array2p(H); p_BT = array2p(BT); p_mass = array2p(mass)
    p_freq = array2p(freq); p_LT = array2p(LT); p_LinvT = array2p(LinvT)
    FL.geometrytransformation_mp_wilsongfmethod_\
        (p_H, p_BT, p_mass, p_freq, p_LT, p_LinvT, byref(c_int(intdim)), byref(c_int(NAtoms)))
    p2array(p_freq, freq); p2array(p_LT, LT); p2array(p_LinvT, LinvT)

# Due to row- and column-major difference, python
#     throws: L^T, B^T
#     fetchs: cartmode^T
def InternalMode2CartesianMode(LT:numpy.ndarray, BT:numpy.ndarray, cartmodeT:numpy.ndarray, intdim:int, cartdim:int) -> None:
    p_LT = array2p(LT); p_BT = array2p(BT)
    p_cartmodeT = array2p(cartmodeT)
    FL.geometrytransformation_mp_internalmode2cartesianmode_\
        (p_LT, p_BT, p_cartmodeT, byref(c_int(intdim)), byref(c_int(cartdim)))
    p2array(p_cartmodeT, cartmodeT)

# ------------------- End --------------------

""" FortranLibrary """
TestFortranLibrary = FL.fortranlibrary_mp_testfortranlibrary_

''' A test program on the python interface of Fortran-Library '''
if __name__ == "__main__":
    print('>>> Testing calling from python... >>>')
    print('\nTime display')
    ShowTime()
    print('\nScientific notation')
    x=3564.1212587; i=0
    x,i=dScientificNotation(x,i)
    print(x - 3.5641212587, i - 3)
    print('\n<<< Calling from python test passed <<<')
    
    print('\n>>> Testing calling from fortran... >>>')
    TestFortranLibrary()
    print('\n<<< Calling from fortran test passed <<<')