#include "ZielonkaTree.hh"
#include <string>

int main(int argc, char** argv) {
    if (argc < 2) return 1;
    ZielonkaTree z("test_condition.txt", std::stoi(argv[1]));
    return 0;
}
