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

from setuptools import setup, find_packages
from setuptools.extension import Extension

# get path to and name of package
name = 'Nsubjettiness'
lname = name.lower()

# function to query a config binary and get the result
fastjet_config = os.environ.get('FASTJET_CONFIG', 'fastjet-config')
def query_config(query):
    return subprocess.check_output([fastjet_config, query]).decode('utf-8').strip()

# get fastjet info
fj_prefix = query_config('--prefix')
fj_cxxflags = query_config('--cxxflags')
fj_ldflags = query_config('--libs')

# get contrib README
with open('README', 'r') as f:
    readme = f.read()

# get contrib version
with open('VERSION', 'r') as f:
    __version__ = f.read().strip() + 'a0'

HELP_MESSAGE = """{name} FastJet Contrib Python Package

Usage: python3 setup.py [COMMAND] [OPTIONS]

Valid options for COMMAND include:
  help - Show this message
  swig - Run SWIG to generate new {lname}.py and Py{name}.cc files; OPTIONS are passed to SWIG
  build_ext - Build the Python extension
  install - Install the Python extension to a standard location
  clean - Remove Python build directories

OPTIONS are passed along to setuptools commands such as build_ext, install, and clean.

Relevant environment variables include:
  FASTJET_CONFIG - Path to fastjet-config binary [defaults to looking for fastjet-config in PATH]
  CXXFLAGS - Compiler flags passed along when building the Python extension module
"""

def show_help():
    print(HELP_MESSAGE.format(name=name, lname=lname))

def run_swig():

    contrib = {'docstring': '{}'.format(readme.replace('"', r'\"')),
               'version': "{}".format(__version__)}

    interface_file = '{}.i'.format(lname)
    template_file = '{}.i.template'.format(lname)
    print('Constructing SWIG interface file {} from {}'.format(interface_file, template_file))

    # read interface template and write interface file
    with open(template_file, 'r') as f_template, open(interface_file, 'w') as f_interface:
        f_interface.write(f_template.read().format(**contrib))

    # form swig options
    fj_swig_interface = os.path.join(fj_prefix, 'share', 'fastjet', 'pyinterface', 'fastjet.i')
    opts = '-fastproxy {} -DFASTJET_SWIG_INTERFACE={}'.format(fj_cxxflags, fj_swig_interface)

    if len(sys.argv) > 2:
        opts += ' ' + ' '.join(sys.argv[2:])

    command = 'swig -python -c++ {} -o Py{}.cc {}.i'.format(opts, name, lname)
    print(command)
    subprocess.run(command.split())

    # move nsubjettiness.py into subdirectory
    os.rename(lname + '.py', os.path.join(lname, lname + '.py'))

def run_setup():

    # get cxxflags from environment, add fastjet cxxflags, and SWIG type table info
    cxxflags = os.environ.get('CXXFLAGS', '').split() + fj_cxxflags.split() + ['-DSWIG_TYPE_TABLE=fastjet']
    setup_path = os.path.abspath(os.path.dirname(__file__))
    ldflags = ['-Wl,-rpath,{}'.format(setup_path)]

    # determine library paths and names for Python
    fj_libdirs, libs = [setup_path], [name]
    for x in fj_ldflags.split():
        if x.startswith('-L'):
            fj_libdirs.append(x[2:])
        elif x.startswith('-l'):
            libs.append(x[2:])
        else:
            ldflags.append(x)

    module = Extension('nsubjettiness._' + lname,
                       sources=['Py{}.cc'.format(name)],
                       language='c++',
                       library_dirs=fj_libdirs,
                       libraries=libs,
                       extra_compile_args=cxxflags,
                       extra_link_args=ldflags)

    setup(
        version=__version__,
        packages=find_packages(),
        ext_modules=[module],
    )

def main():
    commands = {'help': show_help, 'swig': run_swig}
    commands.get(sys.argv[1], run_setup)()

if __name__ == '__main__':
    main()
