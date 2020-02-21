from .basic import *

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
def WilsonBMatrixAndInternalCoordinateq(r:numpy.ndarray, BT:numpy.ndarray, q:numpy.ndarray,\
    cartdim:int, intdim:int) -> None:
    p_r = array2p(r)
    p_BT = array2p(BT); p_q = array2p(q)
    FL.geometrytransformation_mp_wilsonbmatrixandinternalcoordinateq_\
        (p_r, p_BT, p_q, byref(c_int(cartdim)), byref(c_int(intdim)))
    p2array(p_BT, BT); p2array(p_q, q)

# =================== End ===================

# ========== Cartesian <- Internal ==========

def CartesianCoordinater(q:numpy.ndarray, r:numpy.ndarray, intdim:int, cartdim:int,\
    uniquify='none', mass=numpy.array([numpy.nan]), r0=numpy.array([numpy.nan])) -> None:
    p_q = array2p(q)
    p_r = array2p(r)
    n = len(uniquify)
    f = (c_char*n)(); f.value = uniquify.encode('ascii')
    if numpy.isnan(mass[0]): mass=numpy.random.rand(int(cartdim/3))
    if numpy.isnan(r0[0]): r0=numpy.random.rand(cartdim)
    p_mass = array2p(mass); p_r0 = array2p(r0)
    FL.geometrytransformation_mp_cartesiancoordinater_\
        (p_q, p_r, byref(c_int(intdim)), byref(c_int(cartdim)), f, p_mass, p_r0, n)
    p2array(p_r, r)

# =================== End ===================

# ------------------- End --------------------

# --------------- Normal mode ----------------
# Due to row- and column-major difference, python
#     throws: H^T (H^T = H), B^T
#     fetchs: intmode^T, (L^-1)^T, cartmode^T
# Output: freq, intmode^T, (L^-1)^T, cartmode^T
def WilsonGFMethod(H:numpy.ndarray, BT:numpy.ndarray, mass:numpy.ndarray,\
    freq:numpy.ndarray, intmodeT:numpy.ndarray, LinvT:numpy.ndarray, cartmodeT:numpy.ndarray,\
    intdim:int, NAtoms:int) -> None:
    p_H = array2p(H); p_BT = array2p(BT); p_mass = array2p(mass)
    p_freq = array2p(freq); p_intmodeT = array2p(intmodeT); p_LinvT = array2p(LinvT); p_cartmodeT = array2p(cartmodeT)
    FL.geometrytransformation_mp_wilsongfmethod_\
        (p_H, p_BT, p_mass, p_freq, p_intmodeT, p_LinvT, p_cartmodeT,\
        byref(c_int(intdim)), byref(c_int(NAtoms)))
    p2array(p_freq, freq); p2array(p_intmodeT, intmodeT); p2array(p_LinvT, LinvT); p2array(p_cartmodeT, cartmodeT)

# ------------------- End --------------------