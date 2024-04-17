#include "ZielonkaTree.hh"
#include <string>

int main(int argc, char** argv) {
    if (argc < 2) return 1;
    ZielonkaTree z(std::stoi(argv[1]));
    return 0;
}
