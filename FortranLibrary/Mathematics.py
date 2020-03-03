from .basic import *

# No need to fetch constant, directly provide them in python
AMUInAU  = 1822.888486192
AInAU    = 1.8897261339212517
cm_1InAu = 0.000004556335830019422

# -------------------- Quaternion --------------------

def Rotate(q:numpy.ndarray, r:numpy.ndarray) -> None:
    p_q = array2p(q); p_r = array2p(r)
    if len(r.shape) == 1:
        NAtoms = c_int(int(r.shape[0]/3))
    else:
        NAtoms = c_int(r.shape[0])
    FL.mathematics_mp_rotate_(p_q, p_r, byref(NAtoms))
    p2array(p_r, r)

# ----------------------- End ------------------------