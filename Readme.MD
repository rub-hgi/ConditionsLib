# Introduction

This repository contains a C++ library for finding affine decision trees for vectorial Boolean functions
and is published alongside the 2021 Asiacrypt paper "Generic Framework for Key-Guessing Improvements".

# Dependencies

We tested the library with the following dependencies on a Linux machine:

| Dependency                           | Version tested |
|--------------------------------------|----------------|
| [gcc] (https://gcc.gnu.org)          | 11.1.0         |
| [cmake] (https://cmake.org)          | 3.21.2         |
| [NTL] (https://libntl.org)           | 11.5.2         |
| [GMP] (https://gmplib.org)           | 6.2.1          |
| [openmp*] (https://openmp.llvm.org/) | 12.0.1         | 
| [doxygen**] (https://www.doxygen.nl) | 1.9.2          | 

where the dependencies marked with a * are only necessary for compiling all of the examples,
those with ** are only necessary for building the documentation.
We make use of gcc's C++20-support (mainly for std::popcount).

# Building the Library 

The library needs to be compiled (i.e. it is not header only).
It suffices to run cmake on the top-level CMakeLists.txt and then use the resulting makefile for building the library.
```
cd <build directory>
cmake <project's root director>
make all
```

You can install the library system-wide by running `make install` with the appropriate privileges. The exported
target for CMake is `ConditionsLib::ConditionsLib`.

For an example how to use the library with CMake refer to examples-directory.

# Building the Documentation

A pdf file of the code documentation can be found in the doc/-folder.

To build the documentation, run
```
cd <documentation directory>
cmake <project's root directory>
make doxygen
```

# Equivalence Classes of Balanced Functions and 4-bit Sboxes
The equivalence classes of balanced functions 
up to dimension 6 and 4-bit Sboxes in the FunctionData
directory were given to us by Shahram Rasoolzadeh.