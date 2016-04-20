# project template configuration script for basisproject tool

from __future__ import unicode_literals

# ------------------------------------------------------------------------------
# required project files
required = [
  'CMakeLists.txt',
  'BasisProject.cmake'
]

# ------------------------------------------------------------------------------
# optional project files
options = {
  # readme
  'readme' : {
    'desc' : 'Include/exclude README and LICENSE files.',
    'path' : [ 'README.md', 'LICENSE.txt' ]
  }
  # additional configuration files
  'find-mirtk' : {
    'desc' : 'Include/exclude FindMIRTK.cmake file.',
    'path' : [ 'config/FindMIRTK.cmake' ]
  },
  'config-settings' : {
    'desc' : 'Include/exclude custom Settings.cmake file.',
    'path' : [ 'config/Settings.cmake' ]
  },
  'config' : {
    'desc' : 'Include/exclude all custom configuration files.',
    'deps' : [
               'config-settings',
               'config/config.h.in'
             ]
  },
  # source files
  'include' : {
    'desc' : 'Add/remove directory for public header files.',
    'path' : [ 'include/' ]
  },
  'src' : {
    'desc' : 'Add/remove directory for project source files.',
    'path' : [ 'src/CMakeLists.txt' ]
  },
  'tools' : {
    'desc' : 'Add/remove directory for package commands.',
    'path' : [ 'tools/CMakeLists.txt' ]
  },
  # testing tree
  'test' : {
    'desc' : 'Add/remove support for testing.',
    'path' : [
               'test/CMakeLists.txt',
               'test/testClassName.cc'
             ]
  }
}

# ------------------------------------------------------------------------------
# preset template options
presets = {
  'minimal' : {
    'desc' : 'Choose minimal project template.',
    'args' : [ 'include', 'src' ]
  },
  'module' : {
    'desc' : 'Choose internal module template.',
    'args' : [ 'config', 'include', 'src', 'test' ]
  }
  'package' : {
    'desc' : 'Choose external package template.',
    'args' : [ 'readme', 'find-mirtk', 'config', 'include', 'src', 'test' ]
  }
}

# ------------------------------------------------------------------------------
# additional substitutions besides <project>, <template>,...
from datetime import datetime as date
from calendar import month_name, month_abbr
todays = date.today()

substitutions = {
  # fixed computed substitutions
  'date'       : todays.strftime('%x'),
  'day'        : todays.day,
  'month'      : todays.month,
  'month-name' : month_name[todays.month],
  'month-abbr' : month_abbr[todays.month],
  'year'       : todays.year,
  # substitutions which can be overridden using a command option
  'copyright' : {
    'help'    : "Copyright statement optionally including years, but not \"Copyright (c) \" or \". All rights reserved.\".",
    'default' : str(todays.year) + " <author>"
  },
  'contact' : {
    'help'    : "Contact details of person responsible for this MIRTK module.",
    'default' : "<author> <<author>@gmail.com>"
  },
  'license' : {
    'help'    : "Software license.",
    'default' : "Apache License Version 2.0"
  }
}
