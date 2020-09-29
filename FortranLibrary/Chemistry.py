from .basic import *

# Due to row- and column-major difference, python throws mode^T
def Avogadro_Vibration(symbol:List, r:numpy.ndarray, freq:numpy.ndarray, modeT:numpy.ndarray, \
file=Path('avogadro.log')) -> None:
    NAtoms = c_int(len(symbol)); vibdim = c_int(freq.shape[0])
    p_symbol = cast(((c_char*2)*NAtoms.value)(), POINTER(c_char*2))
    for i in range(NAtoms.value):
        p_symbol[i][0] = (symbol[i][0]).encode('ascii')
        if len(symbol[i]) > 1:
            p_symbol[i][1] = symbol[i][1].encode('ascii')
        else:
            p_symbol[i][1] = ' '.encode('ascii')
    p_r = array2p(r); p_freq = array2p(freq); p_modeT = array2p(modeT)
    file_str = str(file); n = len(file_str)
    f = (c_char*n)(); f.value = file_str.encode('ascii')
    func = None
    try:
        func = FL.chemistry_mp_avogadro_vibration_
    except AttributeError:
        func = FL.__chemistry_MOD_avogadro_vibration
    func(byref(NAtoms), p_symbol, p_r, byref(vibdim), p_freq, p_modeT, f, n)

def ghOrthogonalization(grad1:numpy.ndarray, grad2:numpy.ndarray, h:numpy.ndarray) -> None:
    p_grad1 = array2p(grad1); p_grad2 = array2p(grad2); p_h = array2p(h)
    dim = c_int(grad1.shape[0])
    func = None
    try:
        func = FL.chemistry_mp_py_ghorthogonalization_
    except AttributeError:
        func = FL.__chemistry_MOD_py_ghorthogonalization
    func(p_grad1, p_grad2, p_h, byref(dim))
    p2array(p_grad1, grad1); p2array(p_grad2, grad2); p2array(p_h, h)