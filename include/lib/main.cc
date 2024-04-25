#include "ZielonkaTree.hh"
#include <string>

int main(int argc, char** argv) {
    if (argc < 3) return 1;
    ZielonkaTree z(argv[1], std::stoi(argv[2]));
    return 0;
}
