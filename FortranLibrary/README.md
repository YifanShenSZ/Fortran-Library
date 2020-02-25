# Python interface for Fortran-Library
Fortran-Library routines can be called in almost the same manner to fortran. Occasional difference does exist due to python-fortran discrepancy, see details below.

Python interface wraps fortran routines, because:
1. Python recognizes fortran functions by their compiled names:
* Fortran function Func contained in module Mod will be renamed as mod_mp_func_
* You may view the compiled names by `nm libFL.so`
2. Some python data types are immutable:
* For example, built-in types such as: int, float, bool, string, tuple
* Such argument value cannot be changed by passing reference
* Example of mutable data types: numpy.ndarray
3. Optional argument cannot be easily passed to fortran:
* In fortran, function func with optional argument optarg can be called by func(optarg=x)
* When interfaced to python, we can no longer selectively pass optarg to fortran

There are also some python convenience to utilize:
1. Python has specified array length with .shape:
* Fortran routine takes array length as an argument
* It is required in multidimensional case
* Although it can be optional for 1 dimensional, I require it considering stability
2. You are using python, so you must be rich in memory and CPU:
* Fortran code is aimed at memory and CPU demanding jobs
* So in many cases fortran code sacrifices convenience, e.g. using the memory of input argument as work space

So the wrappers:
1. Rename fortran functions back to their origin
2. Return immutable output values in python style
3. Take in optional arguments in python style (actually, default them in wrapper then pass all to fortran)
4. No longer take array length as an argument
5. No longer unnecessarily overwrite input arguments

Weird issue:
1. When a fortran variable length string argument receives a long python string, it will keep only leading 10 characters