#!/bin/sh

# Prepare Paando file structure
mkdir include
mkdir lib
mkdir external


# Install boost libraries
wget https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.gz
tar -xf boost_1_81_0.tar.gz
cd boost_1_81_0
./bootstrap.sh --prefix=../ --with-libraries=system,filesystem,regex
./b2
cd ..
rm -rf boost_1_81_0*

# Install rapidXML
wget https://sourceforge.net/projects/rapidxml/files/rapidxml/rapidxml%201.13/rapidxml-1.13.zip/download
unzip rapidxml-1.13.zip
mv rapidxml-1.13 include/
rm rapidxml-1.13.zip

# Install libconfini
git clone https://github.com/madmurphy/libconfini.git
