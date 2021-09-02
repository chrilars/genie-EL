/* a dummy test program to test the cudd installation */

#include "cuddObj.hh"
#include "dddmp.h"
#include <iostream>
int main() {
    Cudd mgr(0, 0);
    BDD x = mgr.bddVar();
}
