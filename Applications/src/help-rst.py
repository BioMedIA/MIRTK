#!/usr/bin/env python

##############################################################################
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2013-2015 Imperial College London
# Copyright 2013-2015 Andreas Schuh
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

"""
Executes named command with the -h option and reformates the output in reStructuredText.
"""

# ==============================================================================
# Imports
# ==============================================================================

import argparse
import os
import re
import subprocess
import sys
import mirtk

# ==============================================================================
# Auxiliaries
# ==============================================================================

# ------------------------------------------------------------------------------
def get_help(command_name):
    """Get help output of command."""
    return mirtk.check_output([command_name, '-h']).splitlines()

# ------------------------------------------------------------------------------
def get_indent(line):
    """Get number of spaces at start of line."""
    indent = 0
    for c in line:
        if   c == ' ':  indent += 1
        elif c == '\t': indent += 2
        else:           break
    return indent

# ------------------------------------------------------------------------------
def strip_empty_lines(lines):
    """Strip empty lines from start and end of array of lines."""
    empty_line = re.compile('\s*$')
    while len(lines) > 0 and empty_line.match(lines[ 0]): del lines[ 0]
    while len(lines) > 0 and empty_line.match(lines[-1]): del lines[-1]
    return lines

# ------------------------------------------------------------------------------
def join_paragraphs(lines):
    """Join paragraphs separated by empty lines."""
    empty_line            = re.compile('\s*$')
    nonindented_line      = re.compile('^\S')
    rst_directive_pattern = re.compile('.. \w+::|.. \[.+\]\s')
    is_rst_directive      = False
    paragraph = ''
    paragraphs = []
    for line in lines:
        if rst_directive_pattern.match(line):
            if paragraph != '':
                paragraphs.append(paragraph)
            paragraph        = line
            is_rst_directive = True
        elif not is_rst_directive and empty_line.match(line):
            if paragraph != '':
                paragraphs.append(paragraph)
            paragraph        = ''
            is_rst_directive = False
        elif is_rst_directive and nonindented_line.match(line):
            if paragraph != '':
                paragraph = re.sub(r'^\s*\n|\n\s*$', '', paragraph)
                paragraphs.append(paragraph)
            paragraph        = line
            is_rst_directive = False
        elif paragraph == '':
            paragraph = line
        else:
            paragraph = '\n'.join([paragraph, line])
    if paragraph != '':
        if is_rst_directive:
            paragraph = re.sub(r'^\s*\n|\n\s*$', '', paragraph)
        paragraphs.append(paragraph)
    return paragraphs

# ------------------------------------------------------------------------------
def join_arguments(lines):
    """Join help lines corresponding to each argument."""
    empty_line  = re.compile('\s*$')
    start       = re.compile('\s\s\s?\s?-?(\S+)')
    indent      = 0
    option_name = ''
    option_help = []
    arguments   = []
    for line in lines:
        m = start.match(line)
        if m:
            if option_name != '':
                arguments.append((option_name, option_help))
            parts  = line.strip().split(' ' * 3, 1)
            indent = get_indent(line)
            option_name = parts[0]
            option_name = re.sub(r'-\[-\]', '-', option_name)
            m = re.match('-\[no\](\S+)', option_name)
            if m: option_name = '-' + m.group(1) + ', -no' + m.group(1)
            m = re.match('<(.+)>$', option_name)
            if m: option_name = m.group(1)
            if len(parts) == 2:
                option_help = [parts[1].strip()]
            else:
                option_help = []
        elif not empty_line.match(line):
            option_help.append(line.strip())
    if option_name != '':
        arguments.append((option_name, option_help))
    return arguments

# ------------------------------------------------------------------------------
def get_synopsis(output):
    """Get synopsis/usage section from command help output."""
    name  = re.compile('mirtk-(\S)')
    start = re.compile('\s*([sS]ynopsis|[uU]sage)\s*:')
    end   = re.compile('(\S.*:|^)$')
    skip  = True
    lines = []
    for line in output:
        if start.match(line): skip = False
        elif not skip and end.match(line): break
        if not skip:
            line = start.sub('', line, count=1)
            line = name.sub(r'mirtk \1', line)
            lines.append(line.strip())
    lines = strip_empty_lines(lines)
    return lines

# ------------------------------------------------------------------------------
def get_description(output):
    """Get description section from command help output."""
    name  = re.compile('mirtk-(\S)')
    start = re.compile('\s*[dD]escription\s*:\s*$')
    end   = re.compile('\S.*[^:]:\s*$') # must not match reST directives such as ".. math::"
    skip  = True
    lines = []
    for line in output:
        if start.match(line): skip = False
        elif not skip and end.match(line): break
        if not skip:
            line = start.sub('', line, count=1)
            line = name.sub(r'mirtk \1', line)
            lines.append(line)
    if len(lines) == 0:
        synopsis_start = re.compile('\s*([sD]ynopsis|[uU]sage):')
        synopsis_line  = False
        synopsis_done  = False
        for line in output:
            if synopsis_start.match(line): synopsis_line = True
            elif synopsis_line and line == '':
                synopsis_line = False
                synopsis_done = True
            elif synopsis_done:
                if end.match(line): break
                line = start.sub('', line, count=1)
                line = name.sub(r'mirtk \1', line)
                lines.append(line)
    if len(lines) == 0:
        return ["No description available."]
    min_indent = 100
    for line in lines:
        if line == '': continue
        min_indent = min(min_indent, get_indent(line))
    if min_indent < 100:
        for i in range(0, len(lines)):
            lines[i] = re.sub(r' ' * min_indent, '', lines[i], count=1)
    paragraphs = join_paragraphs(lines)
    return paragraphs

# ------------------------------------------------------------------------------
def get_brief_description(description):
    """Get brief description from paragraphs of command description."""
    if description:
        return description[0]
    else:
        return 'No description available.'

# ------------------------------------------------------------------------------
def get_arguments(output):
    """Get help of command arguments."""
    start = re.compile('(^|\S.*)([aA]rguments|[oO]ptions)(.*):\s*$')
    end   = re.compile('\S.*:\s*$')
    section_name  = ''
    section_lines = []
    sections = []
    for line in output:
        m = start.match(line)
        if m:
            line = start.sub('', line, count=1)
            if section_name and section_lines:
                section_args = join_arguments(section_lines)
                sections.append((section_name, section_args))
            section_name = (m.group(1) + m.group(2) + m.group(3)).strip()
            if (re.match(r'[oO]ptional\s+[aA]rguments', section_name) or
                re.match(r'[oO]ptions', section_name)):
                section_name = 'Command options'
            elif (re.match(r'[pP]ositional\s+[aA]rguments', section_name) or
                  re.match(r'[rR]equired\s+[aA]rguments', section_name)):
                section_name = 'Arguments'
            section_lines = []
        elif end.match(line):
            section_name = ''
        if section_name:
            section_lines.append(line)
    if section_name and section_lines:
        section_args = join_arguments(section_lines)
        sections.append((section_name, section_args))
    return sections

# ------------------------------------------------------------------------------
def print_arguments(sections):
    """Print command arguments in reST format."""
    list_item = re.compile('\s*(-|[0-9]+\.)\s')
    for section in sections:
        print('\n\n' + section[0])
        print('-' * len(section[0]))
        for arg in section[1]:
            print('\n.. option:: ' + arg[0] + '\n')
            is_list = False
            for line in arg[1]:
                # insert empty lines before and after each list of items
                if list_item.match(line):
                    if not is_list: print('   ')
                    is_list = True
                elif is_list:
                    print('   ')
                print('   ' + line)

# ------------------------------------------------------------------------------
def get_examples(output):
    """Get examples section from command help output."""
    name  = re.compile('\s*mirtk-(\S)')
    start = re.compile('\s*[eE]xamples:')
    end   = re.compile('.*:$')
    skip  = True
    lines = []
    examples = []
    for line in output:
        if start.match(line):
            line = start.sub('', line, count=1)
            skip = False
        elif not skip and end.match(line):
            break
        if not skip:
            if name.match(line):
                lines = strip_empty_lines(lines)
                if lines:
                    examples.append(lines)
                lines = [name.sub(r'mirtk \1', line)]
            else:
                lines.append(line)
    lines = strip_empty_lines(lines)
    if lines: examples.append(lines)
    return examples

# ==============================================================================
# Main
# ==============================================================================

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    preformatted_indent  = 4
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command_name', help='Name of MIRTK command.')
    parser.add_argument('-output', metavar='file', help='Write help to specified output file. (default: STDOUT)')
    parser.add_argument('-output-brief-description', metavar='file', help='Write brief description to specified output file.')
    parser.add_argument('-synopsis', help='Output command synopsis/brief usage.', action='store_true')
    parser.add_argument('-description', help='Output command description.', action='store_true')
    parser.add_argument('-include-description', help='Include description from specified file.')
    parser.add_argument('-arguments', help='Output command arguments.', action='store_true')
    parser.add_argument('-examples', help='Output command examples.', action='store_true')
    parser.add_argument('-noheaders', help='Do not print section headers.', action='store_true')
    parser.add_argument('-generated', help='Add comment noting that file was generated', action='store_true')
    parser.add_argument('-orphan', help='Add :orphan: in first line to silence warnings if document not included', action='store_true')
    args = parser.parse_args()
    if args.output:
        output_path = os.path.abspath(args.output)
        output_dir  = os.path.dirname(output_path)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        sys.stdout = open(output_path, 'w') 
    # get command help output
    command_name = args.command_name
    command_help = get_help(command_name)
    # split command help output into sections
    if not args.synopsis and not args.description and not args.arguments and not args.examples:
        args.synopsis    = True
        args.description = True
        args.arguments   = True
        args.examples    = True
    if args.include_description: args.description = True
    if args.synopsis:
        synopsis = get_synopsis(command_help)
    else:
        synopsis = ''
    if args.description or args.output_brief_description:
        description = get_description(command_help)
    else:
        description = ''
    # write brief description to file
    if args.output_brief_description:
        output_brief_description_dir = os.path.dirname(args.output_brief_description)
        if not os.path.exists(output_brief_description_dir):
            os.makedirs(output_brief_description_dir)
        output_brief_description = open(args.output_brief_description, 'w')
        if args.generated:
            output_brief_description.write('.. Auto-generated by ' + os.path.basename(sys.argv[0]) +
                     ' from "mirtk ' + command_name + ' -h" output\n\n')
        output_brief_description.write('.. |' + command_name + '-brief-description| replace::\n\n')
        brief_description = get_brief_description(description)
        brief_description = re.sub(r':option:`([^`]*)`', ':option:`' + command_name + ' \\1`', brief_description)
        for line in brief_description.split('\n'):
            output_brief_description.write('   ' + line + '\n')
        output_brief_description.close()
    # print help with reST formatting
    if args.orphan:
        print(':orphan:\n')
    if args.generated:
        print('.. Auto-generated by ' + os.path.basename(sys.argv[0]) +
                ' from "mirtk ' + command_name + ' -h" output\n')
    if not args.noheaders:
        print(command_name)
        print(('=' * len(command_name)))
        print('\n.. program:: ' + command_name)
    # synopsis
    if synopsis:
        if not args.noheaders:
            title = 'Synopsis'
            print('\n\n' + title)
            print(('-' * len(title)) + '\n')
        print('::\n')
        print('    ' + '\n    '.join(synopsis))
    # description
    if not args.noheaders and (args.include_description or description):
        title = 'Description'
        print('\n\n' + title)
        print(('-' * len(title)))
    if args.include_description:
        print('\n.. include:: ' + args.include_description + '\n')
    else:
        raw_preformatted_pattern = re.compile('(^|\n)(\t|  )')
        rst_directive_pattern    = re.compile('.. \w+::|.. \[.+\]\s')
        for paragraph in description:
            if not args.arguments:
                paragraph = re.sub(r':option:`([^`]*)`', ':option:`' + command_name + ' \\1`', paragraph)
            if not rst_directive_pattern.match(paragraph) and raw_preformatted_pattern.search(paragraph):
                print('\n::\n')
                lines = paragraph.split('\n')
                if raw_preformatted_pattern.match(lines[0]):
                    for i in range(0, len(lines)):
                        lines[i] = raw_preformatted_pattern.sub('', lines[i])
                print((' ' * preformatted_indent) + ('\n' + ' ' * preformatted_indent).join(lines))
            else:
                print('\n' + paragraph)
    # arguments
    if args.arguments:
        print_arguments(get_arguments(command_help))
    # examples
    if args.examples: examples = get_examples(command_help)
    else:             examples = ''
    if not args.noheaders and examples:
        title = 'Examples'
        print('\n\n' + title)
        print(('-' * len(title)))
    if examples:
        for i in range(0, len(examples)):
            title = 'Example ' + str(i+1)
            print('\n\n' + title)
            print(('~' * len(title)))
            print('\n**Command**::\n')
            print(' ' * preformatted_indent + examples[i][0])
            if len(examples[i]) > 1:
                print('\n**Output/Description**::\n')
                indent = 100
                for line in examples[i][1:]:
                    indent = min(indent, get_indent(line))
                if indent > 0:
                    for j in range(1, len(examples[i])):
                        examples[i][j] = examples[i][j][indent:]
                print((' ' * preformatted_indent) + ('\n' + ' ' * preformatted_indent).join(examples[i][1:]))
    sys.exit(0)
