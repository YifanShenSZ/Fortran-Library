from .basic import *

ShowTime = FL.general_mp_showtime_

def dScientificNotation(x:float) -> (float, int):
    x_copy = c_double(x); i = c_int(0)
    FL.general_mp_dscientificnotation_(byref(x_copy),byref(i))
    return x_copy.value, i.value