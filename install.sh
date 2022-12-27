# A simple script to install all the programs you need.
#
# Run "./install.sh <path>" to install all programs in "<path>".
# If there is no argument, the <path> is set to "/usr/local/" (maybe sudo will be needed).
#
# Script create a folder "submodules", clone every necessary
# project and install everything in <path>.

# Setting installation directory.
if [ "$1" ];
then
    INSTALLATION_PATH=$1
else
    INSTALLATION_PATH=/usr/local/
fi

rm -r -f submodules # we clean up any existing directory with the same name
mkdir submodules
cd submodules

# installing CUDD
git clone https://github.com/ivmai/cudd
cd cudd
autoreconf -f -i
./configure --enable-shared --enable-obj --enable-dddmp --prefix="$INSTALLATION_PATH"
make
make check
make install
# Now we copy and paste two files
cp config.h $INSTALLATION_PATH/include
cp util/util.h $INSTALLATION_PATH/include
cd ..

# installing Sylvan
git clone --branch v1.6.0 https://github.com/trolando/sylvan
cd sylvan
cmake . -DCMAKE_INSTALL_PREFIX="$INSTALLATION_PATH"
make
make install
cd ..

# cpphoaster
wget 'https://automata.tools/hoa/cpphoafparser/down/cpphoafparser-0.99.2.tgz'
tar xvzf cpphoafparser-0.99.2.tgz
mv cpphoafparser-0.99.2 cpphoafparser
rm cpphoafparser-0.99.2.tgz
cd cpphoafparser
make
# Now we copy and paste folder
cp -r include/* $INSTALLATION_PATH/include
cd ..
