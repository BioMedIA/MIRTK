

## ^^^ Leave two lines blank at top which will be filled by CMake BASIS

##############################################################################
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2013-2016 Imperial College London
# Copyright 2013-2016 Andreas Schuh
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##############################################################################

import os
import sys
import socket
import mirtk

# ============================================================================
# help
# ============================================================================

# ----------------------------------------------------------------------------
def print_help():
    """Print help screen of this executable."""
    print("""
Usage:
  {0} [options] [--] <command> [options of command]
  {0} help <command>...
  {0} [options]

Description:
  This executable is a wrapper for the MIRTK commands. The name of the
  command to execute must be given as first argument.

Arguments:
  command   MIRTK command to execute.

Optional arguments:
  --help, -h      Print help and exit
  --verbose, -v   Print command-line of MIRTK command.
""".format(os.path.basename(__file__)))

# ============================================================================
# main
# ============================================================================

if __name__ == '__main__':
    verbose = 0
    argv    = []
    cmdarg  = False
    if len(sys.argv) < 2:
        print_help()
        sys.exit(1)
    for arg in sys.argv[1:]:
        if cmdarg:
            argv.append(arg)
        elif arg == '--':
            cmdarg = True
        elif arg == '-h' or arg == '-help' or arg == '--help':
            print_help()
            sys.exit(0)
        elif arg == '-v' or arg == '-verbose' or arg == '--verbose':
            verbose += 1
        elif not arg.startswith('-'):
            argv.append(arg)
            cmdarg = True
        else:
            sys.stderr.write('Error: Invalid option: ' + arg)
            sys.exit(1)
    if len(argv) == 0:
        print_help()
        sys.exit(1)
    command = argv[0]
    if command == 'help':
        if len(argv) == 1:
            print_help()
        else:
            for cmd in argv[1:]:
                if mirtk.call([cmd, '-h'], verbose=verbose) != 0:
                    sys.exit(1)
        sys.exit(0)
    else:
        if verbose > 1:
            print('\nHost: ' + socket.gethostname() + '\n')
        exit_code = mirtk.call(argv, verbose=verbose)
        if exit_code != 0:
            sys.stderr.write('Error: ' + command + ' command returned non-zero exit status ' + str(exit_code) + '\n')
        sys.exit(exit_code)
