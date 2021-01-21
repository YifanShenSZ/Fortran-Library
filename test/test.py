'''
This is a test program on Fortran-Library python interface
Correct routines should print close to 0
'''

import FortranLibrary as FL

print(__doc__)

print('\nScientific notation')
x = 3564.1212587; i = 0
x, i = FL.dScientificNotation(x)
print(x - 3.5641212587, i - 3)
