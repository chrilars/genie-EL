
sudo apt-get install cmake
mkdir submodules
cd submodules

# installing CUDD
git clone https://github.com/ivmai/cudd
cd cudd
autoreconf -f -i
./configure --enable-shared --enable-obj --enable-dddmp
make
make check
cd ..

# installing Sylvan
git clone https://github.com/trolando/sylvan
cd sylvan
cmake .
make
cd ..

# cpphoaster
wget 'https://automata.tools/hoa/cpphoafparser/down/cpphoafparser-0.99.2.tgz'
tar xvzf cpphoafparser-0.99.2.tgz
mv cpphoafparser-0.99.2 cpphoafparser
rm cpphoafparser-0.99.2.tgz
cd cpphoafparser
make
cd ..

# todo only for debug version
# googletest
git clone https://github.com/google/googletest
cd googletest
sudo apt-get install libgtest-dev
cmake -Dgtest_build_tests=ON -Dgmock_build_tests=ON .
make
make test
cd ..

cd ..
mkdir build
cd build
cmake ..
make
make install
cd ..

