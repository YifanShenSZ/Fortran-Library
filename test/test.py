from FortranLibrary import *

print('>>> Testing calling from python... >>>')
print('\nTime display')
ShowTime()
print('\nScientific notation')
x=3564.1212587; i=0
x,i=dScientificNotation(x,i)
print(x - 3.5641212587, i - 3)
print('\n<<< Calling from python test passed <<<')

print('\n>>> Testing calling from fortran... >>>')
TestFortranLibrary()
print('\n<<< Calling from fortran test passed <<<')