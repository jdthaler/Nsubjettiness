#!/bin/bash

# some linux-specific CIBW variables
export CIBW_MANYLINUX_X86_64_IMAGE=registry.gitlab.com/pkomiske/fastjet:latest
export CIBW_BUILD="cp*-manylinux_x86_64"

# make contrib shared library
make shared
cp lib*.so /usr/local/lib