# Nsubjettiness Package
#  Questions/Comments?  jthaler@jthaler.net
#    Python questions/comments?  pkomiske@mit.edu
#
#  Copyright (c) 2011-14
#  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
#
#----------------------------------------------------------------------
# This file is part of FastJet contrib.
#
# It is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
#
# It is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this code. If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------

from __future__ import print_function

import os
import subprocess
import sys

from setuptools import setup
from setuptools.extension import Extension

# get path to and name of package
path = os.path.abspath(os.path.dirname(__file__))
name = os.path.basename(path)
lname = name.lower()

# function to query a config binary and get the result
fastjet_config = os.environ.get('FASTJET_CONFIG', 'fastjet-config')
def query_config(query):
    ret = subprocess.run('{} {}'.format(fastjet_config, query).split(), stdout=subprocess.PIPE)
    return str(ret.stdout, 'ascii').strip()

# get fastjet info
fj_prefix = query_config('--prefix')
fj_cxxflags = query_config('--cxxflags')
fj_ldflags = query_config('--libs')

# get contrib README
with open(os.path.join(path, 'README'), 'r') as f:
    readme = f.read()

# get contrib version
with open(os.path.join(path, 'VERSION'), 'r') as f:
    __version__ = f.read().strip()

HELP_MESSAGE = """{name} FastJet Contrib Python Package

Usage: python3 setup.py [COMMAND] [OPTIONS]

Valid options for COMMAND include:
  help - Show this message
  swig - Run SWIG to generate new {name}.py and Py{name}.cc files; OPTIONS are passed to SWIG
  build_ext - Build the Python extension
  install - Install the Python extension to a standard location
  clean - Remove Python build directories

OPTIONS are passed along to setuptools commands such as build_ext, install, and clean.

Relevant environment variables include:
  FASTJET_CONFIG - Path to fastjet-config binary [defaults to looking for fastjet-config in PATH]
  CXXFLAGS - Compiler flags passed along when building the Python extension module
"""

def show_help():
    print(HELP_MESSAGE.format(name=name))

def run_swig():

    contrib = {'docstring': '{}'.format(readme.replace('"', r'\"')),
               'version': "{}".format(__version__)}

    interface_file = '{}.i'.format(lname)
    template_file = '{}.i.template'.format(lname)
    print('Constructing SWIG interface file {} from {}'.format(interface_file, template_file))

    # read interface template and write interface file
    with open(os.path.join(path, template_file), 'r') as f_template, \
         open(os.path.join(path, interface_file), 'w') as f_interface:
        f_interface.write(f_template.read().format(**contrib))

    # form swig options
    fj_swig_interface = os.path.join(fj_prefix, 'share', 'fastjet', 'pyinterface', 'fastjet.i')
    opts = '-fastproxy {} -DSWIG_NUMPY -DFASTJET_SWIG_INTERFACE={}'.format(fj_cxxflags, fj_swig_interface)

    if len(sys.argv) > 2:
        opts += ' ' + ' '.join(sys.argv[2:])

    command = 'swig -python -c++ {} -o Py{}.cc {}.i'.format(opts, name, lname)
    print(command)
    subprocess.run(command.split())

def run_setup():

    # get cxxflags from environment, add fastjet cxxflags, and SWIG type table info
    cxxflags = os.environ.get('CXXFLAGS', '').split() + fj_cxxflags.split() + ['-DSWIG_TYPE_TABLE=fastjet']

    # determine library paths and names for Python
    fj_libdirs = [x[2:] for x in fj_ldflags.split() if x.startswith('-L')]
    libs = [name] + [x[2:] for x in fj_ldflags.split() if x.startswith('-l')]

    module = Extension('_' + lname,
                       sources=['Py{}.cc'.format(name)],
                       language='c++',
                       library_dirs=fj_libdirs,
                       libraries=libs,
                       extra_compile_args=cxxflags)

    setup(
        name=name,
        version=__version__,
        author='Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, TJ Wilkason',
        author_email='jthaler@mit.edu',
        maintainer='Patrick Komiske',
        maintainer_email='pkomiske@mit.edu',
        description='{} FastJet Contrib'.format(name),
        long_description=readme,
        url='https://fastjet.hepforge.org/contrib/',
        py_modules=[lname],
        ext_modules=[module],
        license='GPLv2',
        license_files='COPYING'
    )

def main():
    commands = {'help': show_help, 'swig': run_swig}
    commands.get(sys.argv[1], run_setup)()

if __name__ == '__main__':
    main()
