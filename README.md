# Genie

A library for easily generating experiments on Binary Decision Diagram.


This library uses the following libraries:
- [Cudd](https://github.com/ivmai/cudd)
- [Sylvan](https://github.com/trolando/sylvan)
- [CppHoafParser](https://automata.tools/hoa/cpphoafparser/index.html)

These are necessary to use this library.

# Installation guide

Instructions for easy installation of Genie.

## Installation of Libraries 
In this section we will focus on the installation of the necessary libraries

### Libraries Dependencies
All the necessary standard libraries are listed below.
Sometimes errors caused by missing them are not obvious,
so we advise you to make sure they are there before installation.
- ``git``
- ``g++``
- ``pkg-config``
- ``make``
- ``cmake``
- [Cudd](https://github.com/ivmai/cudd) dependencies:
  - ``autoconf``
  - ``libtool``
- [Sylvan](https://github.com/trolando/sylvan) dependencies:
  - ``libgmp-dev``

###
If you do not have any of: 
[Cudd](https://github.com/ivmai/cudd),
[Sylvan](https://github.com/trolando/sylvan) or 
[CppHoafParser](https://automata.tools/hoa/cpphoafparser/index.html), you can run:
```
./install.sh <path>
```
If you run the script without argument, they set ``<path>`` to ``/usr/local/``.
This script creates a folder ``submodules``, clones every necessary project and installs everything in ``<path>``.
You can also see a ``install.sh`` script to install individual libraries.

## How to run Genie

```
mkdir build
cd build
cmake ..
make
cd ..
```

It is possible that cmake will not be able to find the libraries, then you have to run instead:
```
 cmake .. -DCMAKE_PREFIX_PATH=<semicolon-saperated list of/directories>"
```
where ``<paths>`` is semicolon-seperated list of directories. If you run ``./install.sh <path>`` you can use ``<path>`` instead.

# Generate Documentation

To create documentation run:

```
cd docs
doxygen Doxygen.in
```

and now documentation can be found in `docs/html/index.html`.

# Example of Usage

An example of usage you can find in
[MascotSDS](https://gitlab.mpi-sws.org/kmallik/mascot-sds) or 
[synthesis-with-edge-fairness](https://gitlab.mpi-sws.org/kmallik/synthesis-with-edge-fairness) repository.