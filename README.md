# FairSyn

A library for easily generating experiments on Binary Decision Diagram.

# Install

To run the library we need:
  - [Cudd](https://github.com/ivmai/cudd)
  - [Sylvan](https://github.com/trolando/sylvan)
  - [CppHoafParser](https://automata.tools/hoa/cpphoafparser/index.html)

If you do not have these projects, you can run ``install.sh``. 
This script do not need root access, but you will need to add the appropriate flags to ``cmake``.

If you have installed these projects, you can just run:

```
mkdir build
cd build
cmake ..
make
cd ..
```

It is possible that cmake will not be able to find the libraries, in which case you should add the appropriate flags, for example:
``` cmake .. -DCUDD_HOME="$HOME/FairSyn/submodules/cudd" -DSYLVAN_HOME="$HOME/FairSyn/submodules/sylvan" -DCHP_HOME="$HOME/mpi/cpphoafparser/include" ```

The flags will be saved and the next time you want to run just type ``cmake ..``.