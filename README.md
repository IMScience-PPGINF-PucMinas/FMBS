FMBS
==============

C++ implementation of the fuzzy connection value computation, from the fuzzy marker based segmentation atrticles:

[1] "New hierarchy-based segmentation layer: towards automatic marker proposal." 2021 34th SIBGRAPI Conference on Graphics, Patterns and Images (SIBGRAPI). IEEE, 2021.

[2] "Fuzzy-marker-based segmentation using hierarchies." International Conference on Discrete Geometry and Mathematical Morphology. Cham: Springer International Publishing, 2021.


Installation
------------

**Requires a C++ 14 compiler and cmake**

 - `pip install ./fminmax_cpp`

Build a binary wheel
--------------------
 
A binary wheel ease the redistribution of your project and can be installed with *pip* on a client machine without a compiler.

**Create wheel**

 - `cd fminmax_cpp`
 - `python setup.py bdist_wheel`
 - `pip install ./fminmax_cpp`
 
 The wheel is created in the directory `fminmax_cpp/dist`, it will be named `fminmax-XXXXX.whl` where `XXXXXX` are name tags identifying the current platform and Python version. 
 
**Install wheel**
 
Whells can be installed with *pip*:
 
 - `pip install wheel_name.whl`
 
 Note that a binary wheel is specific to a platform and to a python version (a wheel built on Windows with Python 3.5 can only be installed on Windows with Python 3.5).

Tests
-----

Tests are run automatically at the end of a build: the build will fail if tests are not successful. 

Known Issues
------------

Clang on Linux may not work due to ABI compatibilty issues.
