#! /usr/bin/env python

# ATTENTION: DO NOT use the tokens used by the file update anywhere within
#            this file. Write < basis-custom > instead, for example.

# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2014 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# http://opensource.andreasschuh.com/cmake-basis/download.html#license
# ============================================================================

##
# @file  updatefile.py
# @brief Update file from template file while preserving custom sections (deprecated).
#
# This script is used by the BasisUpdate.cmake module. This module is used to
# update files of a project instantiated from a particular revision of the
# BASIS project template during the configure step of CMake. This way,
# projects pull the changes of the compatible template automatically.
# Sections in the original file which are enclosed by the tokens
# \<basis-custom\> and \</basis-custom\> or \<basis-license\> and
# \</basis-license\> (without trailing spaces) are preserved while all other
# content is replaced by the template file. The customized sections are
# inserted into the template in the order they appear in the original file
# and the template file. If more custom sections are present in the original
# file than in the template file, these custom sections are appended at the
# end of the resulting file.
#
# See the documentation of BasisUpdate.cmake for further details.
#
# @ingroup CMakeHelpers

from __future__ import absolute_import, print_function, unicode_literals

# monkey-patch input to behave like raw_input in Python 2.
# See: http://python3porting.com/differences.html#input-and-raw-input
try:
    input = raw_input
except NameError:
    pass

# modules
import os
import sys
import getopt

# constants
customTag  = "basis-custom"
licenseTag = "basis-license"

tokenCustomStart  = "<"  + customTag  + ">"
tokenCustomEnd    = "</" + customTag  + ">"
tokenLicenseStart = "<"  + licenseTag + ">"
tokenLicenseEnd   = "</" + licenseTag + ">"
tokenKeepTemplate = "REMOVE_THIS_STRING_IF_YOU_WANT_TO_KEEP_YOUR_CHANGES"

# ****************************************************************************
def version (progName):
    """Print version information."""
    print(progName + "1.0.0")

# ****************************************************************************
def usage (progName):
    print("Usage:")
    print("  " + progName + " [options]")
    print()
    print("Required options:")
    print("  [-i --in]       : Filename of original file")
    print("  [-t --template] : Filename of template file")
    print()
    print("Options:")
    print("  [-o --out]     Filename of output file. If this option is not given,")
    print("                 changes are not applied and the exit code can be used")
    print("                 to check whether changes would have been applied.")
    print("  [-f --force]   Force overwrite of output file. Otherwise ask user.")
    print()
    print("Return value:")
    print("  0   Merged output differs from input file and output file was")
    print("      written successfully if option -o or --out was given")
    print("  1   Failed to read or write file")
    print("  2   Nothing changed, input file not overwritten")
    print("  3   Merged output differs from input file but user chose not")
    print("      to overwrite input file")
    print("Example:")
    print("  " + progName + " -i CMakeLists.txt -t CMakeLists.txt.template -o CMakeLists.txt")
    print()
    print("Contact:")
    print("  SBIA Group <sbia-software at uphs.upenn.edu>")

# ****************************************************************************
def help (progName):
    usage (progName)

# ****************************************************************************
def extractCustomizedSections (input):
    custom = []
    start  = 0
    while True:
        start = input.find (tokenCustomStart, start)
        if start == -1:
            break
        start += len (tokenCustomStart)
        next = input.find (tokenCustomStart, start)
        end  = input.find (tokenCustomEnd,   start)
        if end == -1 or (next != -1 and end > next):
            print("WARNING: Found begin of customized section without end token '" + tokenCustomEnd + "'")
            print("WARNING: Will keep template section instead")
            custom.append (tokenKeepTemplate)
        else:
            custom.append (input [start:end])
        start = next
    return custom

# ****************************************************************************
def replaceCustomizedSections (input, custom):
    start  = 0
    end    = -1
    result = input
    for section in custom:
        if start != -1:
            start = result.find (tokenCustomStart, start)
            end   = -1
        if start != -1:
            start += len (tokenCustomStart)
            end = result.find (tokenCustomEnd, start)
        if start == -1 or end == -1:
            if section != tokenKeepTemplate:
                result = result + '\n'
                result = result + '# ' + tokenCustomStart
                result = result + section
                result = result + tokenCustomEnd
        else:
            if section != tokenKeepTemplate:
                result = result [:start] + section + result [end:]
                end    = start + len (section)
        if end != -1:
            start = end + len (tokenCustomEnd)
    return result

# ****************************************************************************
def extractLicenseSections (input):
    license = []
    start   = 0
    while True:
        start = input.find (tokenLicenseStart, start)
        if start == -1:
            break
        start += len (tokenLicenseStart)
        next = input.find (tokenLicenseStart, start)
        end  = input.find (tokenLicenseEnd,   start)
        if end == -1 or (next != -1 and end > next):
            print("WARNING: Found begin of license section without end token '" + tokenLicenseEnd + "'")
            print("WARNING: Will keep template section instead")
            license.append (tokenKeepTemplate)
        else:
            license.append (input [start:end])
        start = next
    return license

# ****************************************************************************
def replaceLicenseSections (input, license):
    start  = 0
    end    = -1
    result = input
    for section in license:
        if start != -1:
            start = result.find (tokenLicenseStart, start)
        if start != -1:
            start += len (tokenLicenseStart)
            end = result.find (tokenLicenseEnd, start)
        if start != -1 and end != -1:
            if section != tokenKeepTemplate:
                result = result [:start] + section + result [end:]
            start = end + len (tokenLicenseEnd)
    return result


# ****************************************************************************
def run (inputFile, templateFile, outputFile, force):
    # open input and output files
    try:
        fIn = open (inputFile, 'r')
    except IOError:
        sys.stderr.write ("Failed to open file '" + inputFile + "'\n")
        return 1
    try:
        fTem = open (templateFile, 'r')
    except IOError:
        sys.stderr.write ("Failed to open file '" + templateFile + "'\n")
        fIn.close ()
        return 1
    # read files
    input = fIn.read ()
    fIn.close ()
    result = fTem.read ()
    fTem.close ()
    # extract custom sections from input and substitute them in template
    sections = extractCustomizedSections (input)
    result   = replaceCustomizedSections (result, sections)
    # extract license sections from input and substitute them in template
    sections = extractLicenseSections (input)
    result   = replaceLicenseSections (result, sections)
    # check whether input file equals output file
    if inputFile == outputFile:
        # check if content changed
        if result == input:
            return 2
        # query user if file should be overwritten
        if not force:
            try:
                sys.stdout.write ("Template of file '" + inputFile + "' has been modified.\n")
                sys.stdout.write ("Do you want to apply the changes (basis-custom sections remain unchanged)? ")
                sys.stdout.flush ()
                apply = input ()
	        if apply != "y" and apply != "yes":
                    return 3
            except EOFError:
                return 1
    # write result to output file if -o option given
    if outputFile != "":
        try:
            fOut = open (outputFile, 'w')
        except IOError:
            sys.stderr.write ("Failed to open file '" + outputFile + "'\n")
            return 1
        fOut.write (result)
        fOut.close ()
    # done
    if result == input:
        return 2
    return 0

# ****************************************************************************
if __name__ == "__main__":
    progName = os.path.basename (sys.argv [0])
    # options
    verbosity    = 0
    inputFile    = ""
    templateFile = ""
    outputFile   = ""
    force        = False
    # get options
    try:
        opts, files = getopt.gnu_getopt (sys.argv [1:], "uhvVfi:t:o:",
          ["usage", "help", "version", "verbose", "force", "in=", "template=", "out="])
    except getopt.GetoptError as err:
        usage (progName)
        print(str(err))
        sys.exit(1)
    # parse command line options
    for o, a in opts:
        if o in ["-V", "--verbose"]:
            verbosity += 1
        elif o in ["-h", "--help", "-u", "--usage"]:
            help (progName)
            sys.exit(0)
        elif o in ["-v", "--version", "--Version"]:
            version (progName)
            sys.exit(0)
        elif o in ["-f", "--force"]:
            force = True
        elif o in ["-i", "--in"]:
            inputFile = os.path.realpath(a)
        elif o in ["-t", "--template"]:
            templateFile = os.path.realpath(a)
        elif o in ["-o", "--out"]:
            outputFile = os.path.realpath(a)
        else:
            assert False, "unhandled option " + o
    # all required inputs specified ?
    if inputFile == "" or templateFile == "":
        usage (progName)
        sys.exit (1)
    # run
    sys.exit (run (inputFile, templateFile, outputFile, force))

