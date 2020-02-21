from .basic import *

def Avogadro_Vibration(NAtoms:int, symbol:numpy.ndarray, r:numpy.ndarray,\
    vibdim:int, freq:numpy.ndarray, mode:numpy.ndarray, FileName='avogadro.log') -> None:
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