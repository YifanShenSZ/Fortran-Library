from .basic import *

ShowTime = FL.general_mp_showtime_

def dScientificNotation(x:float, i:int) -> (float, int):
    x_copy = c_double(x); i_copy = c_int(i)
    FL.general_mp_dscientificnotation_(byref(x_copy),byref(i_copy))
    return x_copy.value, i_copy.value