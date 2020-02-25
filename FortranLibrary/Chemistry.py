from .basic import *

# Due to row- and column-major difference, python throws mode^T
def Avogadro_Vibration(symbol:List, r:numpy.ndarray, freq:numpy.ndarray, modeT:numpy.ndarray, file=Path('avogadro.log')) -> None:
    NAtoms = c_int(len(symbol)); vibdim = c_int(freq.shape[0])
    p_symbol = cast(((c_char*2)*NAtoms)(), POINTER(c_char*2))
    for i in range(NAtoms):
        p_symbol[i][0] = (symbol[i][0]).encode('ascii')
        if len(symbol[i]) < 2:
            p_symbol[i][1] = ' '.encode('ascii')
        else:
            p_symbol[i][1] = symbol[i][1].encode('ascii')
    p_r = array2p(r); p_freq = array2p(freq); p_modeT = array2p(modeT)
    file_str = str(file); n = len(file_str)
    f = (c_char*n)(); f.value = file_str.encode('ascii')
    FL.chemistry_mp_avogadro_vibration_\
        (byref(NAtoms), p_symbol, p_r, byref(vibdim), p_freq, p_modeT, f, n)
