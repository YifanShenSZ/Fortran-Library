# C++ interface for Fortran-Library
Fortran-Library routines can be called in almost the same manner to fortran. Occasional difference does exist due to c++-fortran discrepancy, see details below.

C++ interface wraps fortran routines, because:
1. C++ recognizes fortran functions by their compiled names:
* Fortran function Func contained in module Mod will be renamed as mod_mp_func_ by intel compiler or __mod_MOD_func by GNU compiler
* You may view the compiled names by `nm libFL.so`
2. Optional argument cannot be conveniently passed to fortran:
* In fortran, optional arguments can be trully absent. Function func with optional argument optarg can be called by func(optarg=x) with optarg, or func() without
* In c++, optional arguments are actually all present. Function func with optional argument optarg can be called by func(optarg=x) with user value, or func() with default
* So when interfaced to c++, we can no longer selectively pass optarg to fortran, since c++ passes all arguments

So the wrappers:
1. Rename fortran functions back to their origin
2. Take in optional arguments in c++ style (default them in wrapper then pass all to fortran)

The fortran side also provides some help in 'Interoperability' section:
1. Again, the optional argument issue:
* Some fortran routine's behaviour is controled by the presence of optional arguments
* As noted above, when called from c++, only the 'all present' case will be run
* So other cases now are provided as specific fortran routines

Weird issue:
1. When a fortran variable length string argument receives a long c++ string, it will keep only leading 10 characters