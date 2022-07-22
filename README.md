# FairSyn

A library for easily generating experiments on Binary Decision Diagram.

## Installation

To run the library we need:
  - [Cudd](https://github.com/ivmai/cudd)
  - [Sylvan](https://github.com/trolando/sylvan)
  - [CppHoafParser](https://automata.tools/hoa/cpphoafparser/index.html)

### Libraries Dependencies

- ``c++``
- [Cudd] dependencies:
  - ``pkg-config``
- [Sylvan] dependencies:
  - ``libgmp-dev``

### Installation of Libraries 

If you do not have [Cudd], [Sylvan] or [CppHoafParser], you can run:
```
./install.sh <path>
```
where ``<path>`` is absolute path to folder for installation. 
If you run script without argument, they set ``<path>`` to ``/usr/local/``.
This script create a folder "submodules", clone every necessary project and install everything in ``<path>``.


### How to run FairSyn

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

## Generate Documentation

To create documentation run:

```
cd docs
doxygen Doxygen.in
```

and now documentation can be found in `docs/html/index.html`.