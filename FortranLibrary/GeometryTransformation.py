from .basic import *

class InvolvedMotion:
    def __init__(self, MotionType:str, coeff:float, atom:List):
        self.type  = MotionType
        self.coeff = coeff
        self.atom  = atom.copy()

class IntCoordDef:
    def __init__(self):
        self.motion = []

def StandardizeGeometry(geom:numpy.ndarray, mass:numpy.ndarray,
ref:numpy.ndarray=numpy.array([numpy.nan]), grad:numpy.ndarray=numpy.array([numpy.nan])) -> float:
    p_geom = array2p(geom); p_mass = array2p(mass)
    diff = c_double(0.0)
    NAtoms = c_int(mass.shape[0])
    if len(grad.shape) >= 3:
        NStates = c_int(grad.shape[0])
    else:
        NStates = c_int(1)
    # Determine function name
    func_basic = None
    try:
        func_basic = FL.geometrytransformation_mp_standardizegeometry_basic_
    except AttributeError:
        func_basic = FL.__geometrytransformation_MOD_standardizegeometry_basic
    func_grad = None
    try:
        func_grad = FL.geometrytransformation_mp_standardizegeometry_grad_
    except AttributeError:
        func_grad = FL.__geometrytransformation_MOD_standardizegeometry_grad
    func_ref = None
    try:
        func_ref = FL.geometrytransformation_mp_standardizegeometry_ref_
    except AttributeError:
        func_ref = FL.__geometrytransformation_MOD_standardizegeometry_ref
    func = None
    try:
        func = FL.geometrytransformation_mp_standardizegeometry_
    except AttributeError:
        func = FL.__geometrytransformation_MOD_standardizegeometry
    # Determine which function to run
    if numpy.isnan(ref[0]):
        if numpy.isnan(grad[0]):
            func_basic(p_geom, p_mass, byref(NAtoms))
        else:
            p_grad = array2p(grad)
            func_grad(p_geom, p_mass, byref(NAtoms), byref(NStates), p_grad)
            p2array(p_grad, grad)
    else:
        p_ref = array2p(ref)
        if numpy.isnan(grad[0]):
            func_ref(p_geom, p_mass, byref(NAtoms), p_ref, byref(diff))
        else:
            p_grad = array2p(grad)
            func(p_geom, p_mass, byref(NAtoms), byref(NStates), p_ref, byref(diff), p_grad)
            p2array(p_grad, grad)
    p2array(p_geom, geom)
    return diff.value

# ---------- Cartesian <-> Internal ----------

# We provide a function to load internal coordinate format in python
# instead of fetching GeometryTransformation_IntCoordDef
def FetchInternalCoordinateDefinition(format:str,
file=Path('null')) -> (int, List):
    intdim=0; intcoorddef = []
    if format == 'Columbus7':
        # First line is always 'TEXAS'
        # New internal coordinate line starts with 'K'
        if file.exists():
            with open(file,'r') as f: lines=f.readlines()
        else:
            with open('intcfl','r') as f: lines=f.readlines()
        for i in range(1,len(lines)):
            if lines[i][0] == 'K':
                intdim += 1
                intcoorddef.append(IntCoordDef())
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
            intcoorddef[intdim-1].motion.append(InvolvedMotion(MotionType,coeff,atom))
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
        if file.exists():
            with open(file,'r') as f: lines=f.readlines()
        else:
            with open('IntCoordDef','r') as f: lines=f.readlines()
        for i in range(len(lines)):
            strs = lines[i].split()
            if lines[i][0:6] != '      ':
                intdim += 1
                intcoorddef.append(IntCoordDef())
                strs.pop(0)
            coeff = float(strs[0])
            MotionType = strs[1]
            if MotionType == 'stretching':
                atom = [int(strs[2]), int(strs[3])]
            elif MotionType == 'bending':
                atom = [int(strs[2]), int(strs[3]), int(strs[4])]
            elif MotionType == 'torsion':
                atom = [int(strs[2]), int(strs[3]), int(strs[4]), int(strs[5])]
            elif MotionType == 'OutOfPlane':
                atom = [int(strs[2]), int(strs[3]), int(strs[4]), int(strs[5])]
            intcoorddef[intdim-1].motion.append(InvolvedMotion(MotionType,coeff,atom))
    # Normalized linear combination coefficient
    for i in range(intdim):
        norm2 = 0.0
        for j in range(len(intcoorddef[i].motion)): norm2 += intcoorddef[i].motion[j].coeff * intcoorddef[i].motion[j].coeff
        norm2 = numpy.sqrt(norm2)
        for j in range(len(intcoorddef[i].motion)): intcoorddef[i].motion[j].coeff /= norm2
    return intdim, intcoorddef

def DefineInternalCoordinate(format:str,
file=Path('null')) -> (int, int):
    n1 = len(format)
    f1 = (c_char*n1)(); f1.value = format.encode('ascii')
    if file.exists():
        file_str = str(file)
    else:
        if format == 'Columbus7':
            file_str = 'intcfl'
        else:
            file_str = 'IntCoordDef'
    n2 = len(file_str)
    f2 = (c_char*n2)(); f2.value = file_str.encode('ascii')
    intdim = c_int(0); ID = c_int(0)
    func = None
    try:
        func = FL.geometrytransformation_mp_defineinternalcoordinate_
    except AttributeError:
        func = FL.__geometrytransformation_MOD_defineinternalcoordinate
    func(byref(intdim), byref(ID), f1, f2, n1, n2)
    return intdim.value, ID.value

# ========== Cartesian -> Internal ==========

def InternalCoordinate(r:numpy.ndarray, q:numpy.ndarray,
ID:int=1) -> None:
    p_r = array2p(r)
    p_q = array2p(q)
    cartdim = c_int(r.shape[0]); intdim = c_int(q.shape[0])
    func = None
    try:
        func = FL.geometrytransformation_mp_dinternalcoordinate_
    except AttributeError:
        func = FL.__geometrytransformation_MOD_dinternalcoordinate
    func(p_r, p_q, byref(cartdim), byref(intdim), byref(c_int(ID)))
    p2array(p_q, q)

# Due to row- and column-major difference, python
#     throws: cartgrad^T
#     fetchs:  intgrad^T
def Cartesian2Internal(r:numpy.ndarray, cartgradT:numpy.ndarray, q:numpy.ndarray, intgradT:numpy.ndarray,
ID:int=1) -> None:
    p_r = array2p(r); p_cartgradT = array2p(cartgradT)
    p_q = array2p(q); p_intgradT  = array2p( intgradT)
    cartdim = c_int(r.shape[0]); intdim = c_int(q.shape[0])
    if len(intgradT.shape) >= 3:
        NStates = c_int(intgradT.shape[0])
    else:
        NStates = c_int(1)
    func = None
    try:
        func = FL.geometrytransformation_mp_dcartesian2internal_
    except AttributeError:
        func = FL.__geometrytransformation_MOD_dcartesian2internal
    func(p_r, p_cartgradT, p_q, p_intgradT, byref(cartdim), byref(intdim), byref(NStates), byref(c_int(ID)))
    p2array(p_q, q); p2array(p_intgradT, intgradT)

# Due to row- and column-major difference, python fetchs B^T
def WilsonBMatrixAndInternalCoordinate(r:numpy.ndarray, BT:numpy.ndarray, q:numpy.ndarray,
ID:int=1) -> None:
    p_r = array2p(r)
    p_BT = array2p(BT); p_q = array2p(q)
    cartdim = c_int(r.shape[0]); intdim = c_int(q.shape[0])
    func = None
    try:
        func = FL.geometrytransformation_mp_dwilsonbmatrixandinternalcoordinate_
    except AttributeError:
        func = FL.__geometrytransformation_MOD_dwilsonbmatrixandinternalcoordinate
    func(p_r, p_BT, p_q, byref(cartdim), byref(intdim), byref(c_int(ID)))
    p2array(p_BT, BT); p2array(p_q, q)

# =================== End ===================

# ========== Cartesian <- Internal ==========

def CartesianCoordinate(q:numpy.ndarray, r:numpy.ndarray, \
r0=numpy.array([numpy.nan]), ID:int=1) -> None:
    p_q = array2p(q)
    p_r = array2p(r)
    if numpy.isnan(r0[0]): r0=numpy.random.rand(r.shape[0])
    p_r0 = array2p(r0)
    intdim = c_int(q.shape[0]); cartdim = c_int(r.shape[0])
    func = None
    try:
        func = FL.geometrytransformation_mp_cartesiancoordinate_
    except AttributeError:
        func = FL.__geometrytransformation_MOD_cartesiancoordinate
    func(p_q, p_r, byref(intdim), byref(cartdim), p_r0, byref(c_int(ID)))
    p2array(p_r, r)

# Due to row- and column-major difference, python
#     throws:  intgrad^T
#     fetchs: cartgrad^T
def Internal2Cartesian(q:numpy.ndarray, intgradT:numpy.ndarray, r:numpy.ndarray, cartgradT:numpy.ndarray, \
r0=numpy.array([numpy.nan]), ID:int=1) -> None:
    p_q = array2p(q); p_intgradT  = array2p( intgradT)
    p_r = array2p(r); p_cartgradT = array2p(cartgradT)
    intdim = c_int(q.shape[0]); cartdim = c_int(r.shape[0])
    if len(intgradT.shape) >= 3:
        NStates = c_int(intgradT.shape[0])
    else:
        NStates = c_int(1)
    if numpy.isnan(r0[0]): r0=numpy.random.rand(r.shape[0])
    p_r0 = array2p(r0)
    func = None
    try:
        func = FL.geometrytransformation_mp_internal2cartesian_
    except AttributeError:
        func = FL.__geometrytransformation_MOD_internal2cartesian
    func(p_q, p_intgradT, p_r, p_cartgradT, byref(intdim), byref(cartdim), byref(NStates), p_r0, byref(c_int(ID)))
    p2array(p_r, r); p2array(p_cartgradT, cartgradT)

# =================== End ===================

# ------------------- End --------------------

# --------------- Normal mode ----------------
# Due to row- and column-major difference, python
#     throws: H^T (H^T = H), B^T
#     fetchs: intmode^T, (L^-1)^T, cartmode^T
def WilsonGFMethod(H:numpy.ndarray, BT:numpy.ndarray, mass:numpy.ndarray, \
freq:numpy.ndarray, intmodeT:numpy.ndarray, LinvT:numpy.ndarray, cartmodeT:numpy.ndarray) -> None:
    p_H = array2p(H); p_BT = array2p(BT); p_mass = array2p(mass)
    p_freq = array2p(freq); p_intmodeT = array2p(intmodeT); p_LinvT = array2p(LinvT); p_cartmodeT = array2p(cartmodeT)
    intdim = H.shape[0]; NAtoms = mass.shape[0]
    func = None
    try:
        func = FL.geometrytransformation_mp_wilsongfmethod_
    except AttributeError:
        func = FL.__geometrytransformation_MOD_wilsongfmethod
    func(p_H, p_BT, p_mass, p_freq, p_intmodeT, p_LinvT, p_cartmodeT,\
        byref(c_int(intdim)), byref(c_int(NAtoms)))
    p2array(p_freq, freq); p2array(p_intmodeT, intmodeT); p2array(p_LinvT, LinvT); p2array(p_cartmodeT, cartmodeT)

# ------------------- End --------------------
