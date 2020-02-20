"""
A test program on the python interface of Fortran-Library
"""

import FortranLibrary as FL

print('>>> Testing calling from python... >>>')
print('\nTime display')
FL.ShowTime()
print('\nScientific notation')
x = 3564.1212587; i = 0
x, i = FL.dScientificNotation(x, i)
print(x - 3.5641212587, i - 3)
print('\n<<< Calling from python test passed <<<')

print('\n>>> Testing calling from fortran... >>>')
FL.TestFortranLibrary()
print('\n<<< Calling from fortran test passed <<<')