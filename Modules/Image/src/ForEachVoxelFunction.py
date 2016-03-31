#! /usr/bin/env python

from __future__ import absolute_import, print_function, unicode_literals

import sys

# ------------------------------------------------------------------------------
def to_nary_string(arity):
  return str(arity) + '-ary'
 
# ------------------------------------------------------------------------------
def string_to_arity(arity):
  arity = arity.lower()
  if   arity == 'unary':      return 1
  elif arity == 'binary':     return 2
  elif arity == 'ternary':    return 3
  elif arity == 'quaternary': return 4
  elif arity == 'quinary':    return 5
  elif arity == 'senary':     return 6
  elif arity == 'septenary':  return 7
  elif arity == 'octary':     return 8
  elif arity == 'nonary':     return 9
  return -1

# ------------------------------------------------------------------------------
def arity_to_string(arity):
  if   arity == 1: return 'unary'
  elif arity == 2: return 'binary'
  elif arity == 3: return 'ternary'
  elif arity == 4: return 'quaternary'
  elif arity == 5: return 'quinary'
  elif arity == 6: return 'senary'
  elif arity == 7: return 'septenary'
  elif arity == 8: return 'octary'
  elif arity == 9: return 'nonary'
  return to_nary_string(arity)

# ------------------------------------------------------------------------------
def to_valid_symbol_name(s):
  return s.replace('-', '')
 
# ------------------------------------------------------------------------------
def get_source_name(arity):
  return 'mirtkForEach' + to_valid_symbol_name(arity_to_string(arity)).title() + 'VoxelFunction'

# ------------------------------------------------------------------------------
# parse command-line arguments
if len(sys.argv) != 3:
  print("usage: " + sys.argv[0] + ' <arity> <file>')
  sys.exit(1)
try:
  arity = int(sys.argv[1])
except ValueError:
  arity = string_to_arity(sys.argv[1])
if arity < 1:
  sys.stderr.write('Input argument must be either arity as positive number or a string such as "unary" or "binary"!\n')
  sys.exit(1)
f = open(sys.argv[2], 'w')
if not f:
  sys.stderr.write('Failed to open file ' + sys.argv[2] + '!\n')

# ------------------------------------------------------------------------------
# source file header
source_name = get_source_name(arity)
f.write("""/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 * ATTENTION: This source file has been automatically generated using the code
 *            generator mirtkForEachVoxelFunction.py! This generator is
 *            invoked during CMake configuration of the build system when this
 *            source file is missing from the project.
 *
 *            DO NOT modify this file manually. Instead, modify the code
 *            generator, remove any existing mirtkForEach*VoxelFunction.h
 *            header file from the include/ directory and then re-run CMake.
 *            This will invoke the code generator to re-generate the source files.
 */

""")
# include guard
f.write('#ifndef MIRTK_' + source_name[5:] + '_H\n')
f.write('#define MIRTK_' + source_name[5:] + '_H\n')
# include statements
f.write("""
#include "mirtk/Stream.h"
#include "mirtk/VoxelFunction.h"


namespace mirtk {

""")

# ------------------------------------------------------------------------------
# generate ForEach function for each combination of const and non-const image arguments
def gencode(arity, num_const):
  arity_string = arity_to_string(arity).title()
  cfg = {}
  # use last image as reference for image size and inside check as it is usually
  # the output image
  cfg['refim'] = 'im' + str(arity)
  cfg['refp']  = 'p'  + str(arity)
  # other settings
  if   num_const == 0:     cfg['num_const_comment']  = str(arity) + ' non-const'
  elif num_const == arity: cfg['num_const_comment']  = str(arity) + ' const'
  else:                    cfg['num_const_comment']  = str(num_const) + ' const, ' + str(arity - num_const) + ' non-const'
  if arity > 1:            cfg['num_const_comment'] += ' images'
  else:                    cfg['num_const_comment'] += ' image'
  cfg['class_name1']           = arity_string + 'ForEachVoxelBody'
  cfg['class_name2']           = arity_string + 'ForEachVoxelIfBody'
  if num_const > 0:
    cfg['class_name1']        += '_'
    cfg['class_name2']        += '_'
    if num_const < arity:
      cfg['class_name1']        += str(num_const)
      cfg['class_name2']        += str(num_const)
    cfg['class_name1']        += 'Const'
    cfg['class_name2']        += 'Const'
  cfg['member_declaration']    = ''
  cfg['class_T']               = ''
  cfg['T']                     = ''
  cfg['constructor_args1']     = ''
  cfg['constructor_args2']     = ''
  cfg['init_list']             = ''
  cfg['copy_list']             = ''
  cfg['init_pointers']         = ''
  cfg['init_pointers_1D']      = ''
  cfg['init_pointers_2D']      = ''
  cfg['init_pointers_3D']      = ''
  cfg['preincrement_pointers'] = ''
  cfg['inc_pointers_col']      = ''
  cfg['inc_pointers_row']      = ''
  cfg['inc_pointers_page']     = ''
  cfg['pargs']                 = ''
  cfg['imparams_by_reference'] = ''
  cfg['imargs']                = ''
  cfg['impargs']               = ''
  for i in range(1, arity+1):
    n = str(i)
    if i > 1:
      cfg['member_declaration']    += '\n' + (' ' * 2)
      cfg['class_T']               += ', '
      cfg['T']                     += ', '
      cfg['constructor_args1']     += ',\n' + (' ' * (2 + len(cfg['class_name1']) + 1))
      cfg['constructor_args2']     += ',\n' + (' ' * (2 + len(cfg['class_name2']) + 1))
      cfg['init_list']             += ', '
      cfg['copy_list']             += ', '
      cfg['init_pointers']         += '\n' + (' ' * 4)
      cfg['init_pointers_1D']      += '\n' + (' ' * 4)
      cfg['init_pointers_2D']      += '\n' + (' ' * 4)
      cfg['init_pointers_3D']      += '\n' + (' ' * 4)
      cfg['preincrement_pointers'] += ', '
      cfg['inc_pointers_col']      += ', '
      cfg['inc_pointers_row']      += ', '
      cfg['inc_pointers_page']     += ', '
      cfg['pargs']                 += ', '
      cfg['imparams_by_reference'] += ', '
      cfg['imargs']                += ', '
      cfg['impargs']               += ', '
    if num_const > 0:
      if i <= num_const:
        cfg['member_declaration']    += 'const '
        cfg['constructor_args1']     += 'const '
        cfg['constructor_args2']     += 'const '
        cfg['init_pointers']         += 'const '
        cfg['init_pointers_1D']      += 'const '
        cfg['init_pointers_2D']      += 'const '
        cfg['init_pointers_3D']      += 'const '
        cfg['imparams_by_reference'] += 'const '
      else:
        cfg['member_declaration']  += '      '
        cfg['constructor_args1']   += '      '
        cfg['constructor_args2']   += '      '
        cfg['init_pointers']       += '      '
        cfg['init_pointers_1D']    += '      '
        cfg['init_pointers_2D']    += '      '
        cfg['init_pointers_3D']    += '      '
    cfg['member_declaration']      += 'GenericImage<T' + n + '> &im' + n + ';'
    cfg['imparams_by_reference']   += 'GenericImage<T' + n + '> &im' + n
    cfg['class_T']                 += 'class T' + n
    cfg['T']                       +=       'T' + n
    cfg['constructor_args1']       += 'GenericImage<T' + n + '> &im' + n
    cfg['constructor_args2']       += 'GenericImage<T' + n + '> &im' + n
    cfg['init_list']               += 'im' + n + '(im' + n + ')'
    cfg['copy_list']               += 'im' + n + '(o.im' + n + ')'
    cfg['init_pointers']           += 'T' + n + ' *p' + n + ' = im' + n + '.IsEmpty() ? NULL : im' + n + '.GetPointerToVoxels();'
    cfg['init_pointers_1D']        += 'T' + n + ' *p' + n + ' = im' + n + '.IsEmpty() ? NULL : im' + n + '.GetPointerToVoxels() + re.begin();'
    cfg['init_pointers_2D']        += 'T' + n + ' *p' + n + ' = im' + n + '.IsEmpty() ? NULL : im' + n + '.GetPointerToVoxels(bi, bj, this->_k, this->_l);'
    cfg['init_pointers_3D']        += 'T' + n + ' *p' + n + ' = im' + n + '.IsEmpty() ? NULL : im' + n + '.GetPointerToVoxels(bi, bj, bk, this->_l);'
    cfg['preincrement_pointers']   += '++p' + n
    cfg['inc_pointers_col']        += 'p' + n + ' +=  1'
    cfg['inc_pointers_row']        += 'p' + n + ' += s1'
    cfg['inc_pointers_page']       += 'p' + n + ' += s2'
    cfg['pargs']                   += 'p' + n
    cfg['imargs']                  += 'im' + n
    cfg['impargs']                 += '*im' + n
  cfg['constructor_args1']         += ',\n' + (' ' * (2 + len(cfg['class_name1']) + 1))
  cfg['constructor_args2']         += ',\n' + (' ' * (2 + len(cfg['class_name2']) + 1))
  cfg['constructor_args1']         += 'VoxelFunc &vf'
  cfg['constructor_args2']         += 'VoxelFunc &vf, OutsideFunc &of'
  cfg['imparams_by_pointer'] = cfg['imparams_by_reference'].replace('&', '*')
  cfg['assert_is_not_reduction']         = "if (VoxelFunc::IsReduction()) _foreach%svoxelfunction_must_not_be_reduction();"                               % arity_string.lower()
  cfg['assert_neither_is_not_reduction'] = "if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreach%svoxelfunction_must_not_be_reduction();" % arity_string.lower()
  f.write("""
// =============================================================================
// %(num_const_comment)s
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of %(num_const_comment)s
 */
template <%(class_T)s, class VoxelFunc>
struct %(class_name1)s : public ForEachVoxelBody<VoxelFunc>
{
  %(member_declaration)s

  /// Constructor
  %(class_name1)s(%(constructor_args1)s)
  :
    ForEachVoxelBody<VoxelFunc>(vf, im1.Attributes()), %(init_list)s
  {}

  /// Copy constructor
  %(class_name1)s(const %(class_name1)s &o)
  :
    ForEachVoxelBody<VoxelFunc>(o), %(copy_list)s
  {}

  /// Split constructor
  %(class_name1)s(%(class_name1)s &o, split s)
  :
    ForEachVoxelBody<VoxelFunc>(o, s), %(copy_list)s
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    %(init_pointers)s

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, %(preincrement_pointers)s) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<%(class_name1)s *>(this)->_VoxelFunc(i, j, k, l, %(pargs)s);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    %(init_pointers_1D)s

    for (int idx = re.begin(); idx < re.end(); ++idx, %(inc_pointers_col)s) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<%(class_name1)s *>(this)->_VoxelFunc(%(refim)s, idx, %(pargs)s);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = %(refim)s.GetX() - (ei - bi);

    %(init_pointers_2D)s

    for (int j = bj; j < ej; ++j, %(inc_pointers_row)s)
    for (int i = bi; i < ei; ++i, %(inc_pointers_col)s) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<%(class_name1)s *>(this)->_VoxelFunc(i, j, this->_k, this->_l, %(pargs)s);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  %(refim)s.GetX() - (ei - bi);
    const int s2 = (%(refim)s.GetY() - (ej - bj)) * %(refim)s.GetX();

    %(init_pointers_3D)s

    for (int k = bk; k < ek; ++k, %(inc_pointers_page)s)
    for (int j = bj; j < ej; ++j, %(inc_pointers_row)s)
    for (int i = bi; i < ei; ++i, %(inc_pointers_col)s) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<%(class_name1)s *>(this)->_VoxelFunc(i, j, k, this->_l, %(pargs)s);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of %(num_const_comment)s
 */
template <%(class_T)s,
          class VoxelFunc, class OutsideFunc = NaryVoxelFunction::NOP,
          class Domain = ForEachVoxelDomain::Foreground>
struct %(class_name2)s : public ForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  %(member_declaration)s

  /// Constructor
  %(class_name2)s(%(constructor_args2)s)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of, im1.Attributes()), %(init_list)s
  {}

  /// Copy constructor
  %(class_name2)s(const %(class_name2)s &o)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), %(copy_list)s
  {}

  /// Split constructor
  %(class_name2)s(%(class_name2)s &o, split s)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), %(copy_list)s
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    %(init_pointers)s

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, %(preincrement_pointers)s) {
      if (Domain::IsInside(%(refim)s, i, j, k, l, %(refp)s)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<%(class_name2)s *>(this)->_VoxelFunc  (i, j, k, l, %(pargs)s);
      } else const_cast<%(class_name2)s *>(this)->_OutsideFunc(i, j, k, l, %(pargs)s);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    %(init_pointers_1D)s

    for (int idx = re.begin(); idx < re.end(); ++idx, %(inc_pointers_col)s) {
      if (Domain::IsInside(%(refim)s, idx, %(refp)s)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<%(class_name2)s *>(this)->_VoxelFunc  (%(refim)s, idx, %(pargs)s);
      } else const_cast<%(class_name2)s *>(this)->_OutsideFunc(%(refim)s, idx, %(pargs)s);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = %(refim)s.GetX() - (ei - bi);

    %(init_pointers_2D)s

    for (int j = bj; j < ej; ++j, %(inc_pointers_row)s)
    for (int i = bi; i < ei; ++i, %(inc_pointers_col)s) {
      if (Domain::IsInside(%(refim)s, i, j, this->_k, this->_l, %(refp)s)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<%(class_name2)s *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, %(pargs)s);
      } else const_cast<%(class_name2)s *>(this)->_OutsideFunc(i, j, this->_k, this->_l, %(pargs)s);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  %(refim)s.GetX() - (ei - bi);
    const int s2 = (%(refim)s.GetY() - (ej - bj)) * %(refim)s.GetX();

    %(init_pointers_3D)s

    for (int k = bk; k < ek; ++k, %(inc_pointers_page)s)
    for (int j = bj; j < ej; ++j, %(inc_pointers_row)s)
    for (int i = bi; i < ei; ++i, %(inc_pointers_col)s) {
      if (Domain::IsInside(%(refim)s, i, j, k, this->_l, %(refp)s)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<%(class_name2)s *>(this)->_VoxelFunc  (i, j, k, this->_l, %(pargs)s);
      } else const_cast<%(class_name2)s *>(this)->_OutsideFunc(i, j, k, this->_l, %(pargs)s);
    }
  }
};

// -----------------------------------------------------------------------------
// ForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachScalar(%(imparams_by_pointer)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(impargs)s, vf);
  blocked_range<int> re(0, %(refim)s->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ForEachScalar(%(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(%(imparams_by_pointer)s, VoxelFunc &vf)
{
  if (%(refim)s->GetTSize()) {
    ForEachScalar(%(impargs)s, vf);
  } else {
    %(class_name1)s<%(T)s, VoxelFunc> body(%(impargs)s, vf);
    blocked_range<int> re(0, %(refim)s->GetNumberOfVoxels() / %(refim)s->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ForEachVoxel(%(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, %(imparams_by_pointer)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(impargs)s, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ForEachVoxel(attr, %(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, %(imparams_by_pointer)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(impargs)s, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ForEachVoxel(re, %(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, %(imparams_by_pointer)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(impargs)s, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ForEachVoxel(re, %(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, %(imparams_by_pointer)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(impargs)s, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ForEachVoxel(re, %(impargs)s, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachScalar(%(imparams_by_reference)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(imargs)s, vf);
  blocked_range<int> re(0, %(refim)s.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ForEachScalar(%(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(%(imparams_by_reference)s, VoxelFunc &vf)
{
  if (%(refim)s.GetTSize()) {
    ForEachScalar(%(imargs)s, vf);
  } else {
    %(class_name1)s<%(T)s, VoxelFunc> body(%(imargs)s, vf);
    blocked_range<int> re(0, %(refim)s.GetNumberOfVoxels() / %(refim)s.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ForEachVoxel(%(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, %(imparams_by_reference)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(imargs)s, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ForEachVoxel(attr, %(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, %(imparams_by_reference)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(imargs)s, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ForEachVoxel(re, %(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, %(imparams_by_reference)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(imargs)s, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ForEachVoxel(re, %(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, %(imparams_by_reference)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(imargs)s, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ForEachVoxel(re, %(imargs)s, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(%(imparams_by_pointer)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(impargs)s, vf, of);
  blocked_range<int> re(0, %(refim)s->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, %(imparams_by_pointer)s)
{
  %(assert_neither_is_not_reduction)s
  ForEachScalarIf<Domain>(%(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachScalarIf(%(imparams_by_pointer)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(%(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ForEachScalarIf<Domain>(%(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(%(imparams_by_pointer)s, VoxelFunc &vf, OutsideFunc &of)
{
  if (%(refim)s->GetTSize()) {
    ForEachScalarIf<Domain>(%(impargs)s, vf, of);
  } else {
    %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(impargs)s, vf, of);
    blocked_range<int> re(0, %(refim)s->GetNumberOfVoxels() / %(refim)s->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, %(imparams_by_pointer)s)
{
  %(assert_neither_is_not_reduction)s
  ForEachVoxelIf<Domain>(%(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachVoxelIf(%(imparams_by_pointer)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(%(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ForEachVoxelIf<Domain>(%(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, %(imparams_by_pointer)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(impargs)s, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, %(imparams_by_pointer)s)
{
  %(assert_neither_is_not_reduction)s
  ForEachVoxelIf<Domain>(attr, %(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, %(imparams_by_pointer)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, %(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ForEachVoxelIf<Domain>(attr, %(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, %(imparams_by_pointer)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(impargs)s, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, %(imparams_by_pointer)s)
{
  %(assert_neither_is_not_reduction)s
  ForEachVoxelIf<Domain>(re, %(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, %(imparams_by_pointer)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(impargs)s, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, %(imparams_by_pointer)s)
{
  %(assert_neither_is_not_reduction)s
  ForEachVoxelIf<Domain>(re, %(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, %(imparams_by_pointer)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, %(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ForEachVoxelIf<Domain>(re, %(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, %(imparams_by_pointer)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(impargs)s, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, %(imparams_by_pointer)s)
{
  %(assert_neither_is_not_reduction)s
  ForEachVoxelIf<Domain>(%(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, %(imparams_by_pointer)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, %(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ForEachVoxelIf<Domain>(re, %(impargs)s, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(%(imparams_by_reference)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(imargs)s, vf, of);
  blocked_range<int> re(0, %(refim)s.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, %(imparams_by_reference)s)
{
  %(assert_neither_is_not_reduction)s
  ForEachScalarIf<Domain>(%(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachScalarIf(%(imparams_by_reference)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(%(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ForEachScalarIf<Domain>(%(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(%(imparams_by_reference)s, VoxelFunc &vf, OutsideFunc &of)
{
  if (%(refim)s.GetTSize()) {
    ForEachVoxelIf<Domain>(%(imargs)s, vf, of);
  } else {
    %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(imargs)s, vf, of);
    blocked_range<int> re(0, %(refim)s.GetNumberOfVoxels() / %(refim)s.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, %(imparams_by_reference)s)
{
  %(assert_neither_is_not_reduction)s
  ForEachVoxelIf<Domain>(%(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachVoxelIf(%(imparams_by_reference)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(%(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ForEachVoxelIf<Domain>(%(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, %(imparams_by_reference)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(imargs)s, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, %(imparams_by_reference)s)
{
  %(assert_neither_is_not_reduction)s
  ForEachVoxelIf<Domain>(attr, %(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, %(imparams_by_reference)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, %(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ForEachVoxelIf<Domain>(attr, %(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, %(imparams_by_reference)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(imargs)s, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, %(imparams_by_reference)s)
{
  %(assert_neither_is_not_reduction)s
  ForEachVoxelIf<Domain>(re, %(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, %(imparams_by_reference)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(imargs)s, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, %(imparams_by_reference)s)
{
  %(assert_neither_is_not_reduction)s
  ForEachVoxelIf<Domain>(re, %(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, %(imparams_by_reference)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, %(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ForEachVoxelIf<Domain>(re, %(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, %(imparams_by_reference)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(imargs)s, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, %(imparams_by_reference)s)
{
  %(assert_neither_is_not_reduction)s
  ForEachVoxelIf<Domain>(re, %(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, %(imparams_by_reference)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, %(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ForEachVoxelIf<Domain>(re, %(imargs)s, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachScalar(%(imparams_by_pointer)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(impargs)s, vf);
  blocked_range<int> re(0, %(refim)s->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachScalar(%(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(%(imparams_by_pointer)s, VoxelFunc &vf)
{
  if (%(refim)s->GetTSize()) {
    ParallelForEachScalar(%(impargs)s, vf);
  } else {
    %(class_name1)s<%(T)s, VoxelFunc> body(%(impargs)s, vf);
    blocked_range<int> re(0, %(refim)s->GetNumberOfVoxels() / %(refim)s->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxel(%(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, %(imparams_by_pointer)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(impargs)s, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxel(attr, %(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, %(imparams_by_pointer)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(impargs)s, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxel(re, %(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, %(imparams_by_pointer)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(impargs)s, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxel(re, %(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, %(imparams_by_pointer)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(impargs)s, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxel(re, %(impargs)s, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachScalar(%(imparams_by_reference)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(imargs)s, vf);
  blocked_range<int> re(0, %(refim)s.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachScalar(%(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(%(imparams_by_reference)s, VoxelFunc &vf)
{
  if (%(refim)s.GetTSize()) {
    ParallelForEachScalar(%(imargs)s, vf);
  } else {
    %(class_name1)s<%(T)s, VoxelFunc> body(%(imargs)s, vf);
    blocked_range<int> re(0, %(refim)s.GetNumberOfVoxels() / %(refim)s.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxel(%(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, %(imparams_by_reference)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(imargs)s, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxel(attr, %(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, %(imparams_by_reference)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(imargs)s, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxel(re, %(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, %(imparams_by_reference)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(imargs)s, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxel(re, %(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, %(imparams_by_reference)s, VoxelFunc &vf)
{
  %(class_name1)s<%(T)s, VoxelFunc> body(%(imargs)s, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <%(class_T)s, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxel(re, %(imargs)s, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(%(imparams_by_pointer)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(impargs)s, vf, of);
  blocked_range<int> re(0, %(refim)s->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, %(imparams_by_pointer)s)
{
  %(assert_neither_is_not_reduction)s
  ParallelForEachScalarIf<Domain>(%(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachScalarIf(%(imparams_by_pointer)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(%(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachScalarIf<Domain>(%(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(%(imparams_by_pointer)s, VoxelFunc &vf, OutsideFunc &of)
{
  if (%(refim)s->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(%(impargs)s, vf, of);
  } else {
    %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(impargs)s, vf, of);
    blocked_range<int> re(0, %(refim)s->GetNumberOfVoxels() / %(refim)s->GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, %(imparams_by_pointer)s)
{
  %(assert_neither_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(%(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(%(imparams_by_pointer)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(%(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(%(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, %(imparams_by_pointer)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(impargs)s, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, %(imparams_by_pointer)s)
{
  %(assert_neither_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(attr, %(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, %(imparams_by_pointer)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, %(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(attr, %(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, %(imparams_by_pointer)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(impargs)s, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, %(imparams_by_pointer)s)
{
  %(assert_neither_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(re, %(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, %(imparams_by_pointer)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, %(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(re, %(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, %(imparams_by_pointer)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(impargs)s, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, %(imparams_by_pointer)s)
{
  %(assert_neither_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(re, %(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, %(imparams_by_pointer)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, %(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(re, %(impargs)s, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, %(imparams_by_pointer)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(impargs)s, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, %(imparams_by_pointer)s)
{
  %(assert_neither_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(re, %(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, %(imparams_by_pointer)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, %(impargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, %(imparams_by_pointer)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(re, %(impargs)s, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(%(imparams_by_reference)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(imargs)s, vf, of);
  blocked_range<int> re(0, %(refim)s.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, %(imparams_by_reference)s)
{
  %(assert_neither_is_not_reduction)s
  ParallelForEachScalarIf<Domain>(%(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachScalarIf(%(imparams_by_reference)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(%(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachScalarIf<Domain>(%(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(%(imparams_by_reference)s, VoxelFunc &vf, OutsideFunc &of)
{
  if (%(refim)s.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(%(imargs)s, vf, of);
  } else {
    %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(imargs)s, vf, of);
    blocked_range<int> re(0, %(refim)s.GetNumberOfVoxels() / %(refim)s.GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, %(imparams_by_reference)s)
{
  %(assert_neither_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(%(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(%(imparams_by_reference)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(%(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(%(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, %(imparams_by_reference)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(imargs)s, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, %(imparams_by_reference)s)
{
  %(assert_neither_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(attr, %(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, %(imparams_by_reference)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, %(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(attr, %(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, %(imparams_by_reference)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(imargs)s, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, %(imparams_by_reference)s)
{
  %(assert_neither_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(re, %(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, %(imparams_by_reference)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, %(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(re, %(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, %(imparams_by_reference)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(imargs)s, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, %(imparams_by_reference)s)
{
  %(assert_neither_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(re, %(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, %(imparams_by_reference)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, %(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(re, %(imargs)s, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, %(imparams_by_reference)s, VoxelFunc &vf, OutsideFunc &of)
{
  %(class_name2)s<%(T)s, VoxelFunc, OutsideFunc, Domain> body(%(imargs)s, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, %(imparams_by_reference)s)
{
  %(assert_neither_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(re, %(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, %(imparams_by_reference)s, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, %(imargs)s, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, %(class_T)s, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, %(imparams_by_reference)s)
{
  %(assert_is_not_reduction)s
  ParallelForEachVoxelIf<Domain>(re, %(imargs)s, vf);
}
""" % cfg)

# ------------------------------------------------------------------------------
# main
f.write("""
inline void _foreach%svoxelfunction_must_not_be_reduction()
{
  cerr << "(Parallel)ForEachVoxel(If): Voxel reductions must be passed by reference!"
               " Pass voxel functor object(s) as last argument(s) instead of first." << endl;
  exit(1);
}

""" % arity_to_string(arity))

num_const = arity
while num_const >= 0:
  gencode(arity, num_const)
  num_const -= 1

# ------------------------------------------------------------------------------
# footer - end of namespace and include guard
f.write("""

} // namespace mirtk

#endif
""")
f.close()
