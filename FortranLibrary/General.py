from .basic import *

ShowTime = None
try:
    ShowTime = FL.general_mp_showtime_
except AttributeError:
    ShowTime = FL.__general_MOD_showtime

def dScientificNotation(x:float) -> (float, int):
    x_copy = c_double(x); i = c_int(0)
    func = None
    try:
        func = FL.general_mp_dscientificnotation_
    except AttributeError:
        func = FL.__general_MOD_dscientificnotation
    func(byref(x_copy), byref(i))
    return x_copy.value, i.value