#!/bin/bash

# fastjet prereqs
brew upgrade cgal

# install fastjet
git clone https://gitlab.com/pkomiske/fastjet.git
cd fastjet
git submodule init
git submodule update
autoreconf -i
export CXXFLAGS=-std=c++14
./configure --prefix=/usr/local --enable-pyext --enable-cgal --enable-cgal-header-only --disable-monolithic --disable-allplugins --disable-debug PYTHON=python3 PYTHON_CONFIG=python3-config
make -j2 install
cd ..

# make contrib shared library
make shared
cp lib*.dylib /usr/local/lib