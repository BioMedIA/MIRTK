/*
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

#ifndef MIRTK_ForEachQuaternaryVoxelFunction_H
#define MIRTK_ForEachQuaternaryVoxelFunction_H

#include "mirtk/Stream.h"
#include "mirtk/VoxelFunction.h"


namespace mirtk {


inline void _foreachquaternaryvoxelfunction_must_not_be_reduction()
{
  cerr << "(Parallel)ForEachVoxel(If): Voxel reductions must be passed by reference!"
               " Pass voxel functor object(s) as last argument(s) instead of first." << endl;
  exit(1);
}


// =============================================================================
// 4 const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 4 const images
 */
template <class T1, class T2, class T3, class T4, class VoxelFunc>
struct QuaternaryForEachVoxelBody_Const : public ForEachVoxelBody<VoxelFunc>
{
  const GenericImage<T1> &im1;
  const GenericImage<T2> &im2;
  const GenericImage<T3> &im3;
  const GenericImage<T4> &im4;

  /// Constructor
  QuaternaryForEachVoxelBody_Const(const GenericImage<T1> &im1,
                                   const GenericImage<T2> &im2,
                                   const GenericImage<T3> &im3,
                                   const GenericImage<T4> &im4,
                                   VoxelFunc &vf)
  :
    ForEachVoxelBody<VoxelFunc>(vf, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4)
  {}

  /// Copy constructor
  QuaternaryForEachVoxelBody_Const(const QuaternaryForEachVoxelBody_Const &o)
  :
    ForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Split constructor
  QuaternaryForEachVoxelBody_Const(QuaternaryForEachVoxelBody_Const &o, split s)
  :
    ForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody_Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody_Const *>(this)->_VoxelFunc(im4, idx, p1, p2, p3, p4);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im4.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody_Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4);
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

    const int s1 =  im4.GetX() - (ei - bi);
    const int s2 = (im4.GetY() - (ej - bj)) * im4.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody_Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 4 const images
 */
template <class T1, class T2, class T3, class T4,
          class VoxelFunc, class OutsideFunc = NaryVoxelFunction::NOP,
          class Domain = ForEachVoxelDomain::Foreground>
struct QuaternaryForEachVoxelIfBody_Const : public ForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const GenericImage<T1> &im1;
  const GenericImage<T2> &im2;
  const GenericImage<T3> &im3;
  const GenericImage<T4> &im4;

  /// Constructor
  QuaternaryForEachVoxelIfBody_Const(const GenericImage<T1> &im1,
                                     const GenericImage<T2> &im2,
                                     const GenericImage<T3> &im3,
                                     const GenericImage<T4> &im4,
                                     VoxelFunc &vf, OutsideFunc &of)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4)
  {}

  /// Copy constructor
  QuaternaryForEachVoxelIfBody_Const(const QuaternaryForEachVoxelIfBody_Const &o)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Split constructor
  QuaternaryForEachVoxelIfBody_Const(QuaternaryForEachVoxelIfBody_Const &o, split s)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4) {
      if (Domain::IsInside(im4, i, j, k, l, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      if (Domain::IsInside(im4, idx, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (im4, idx, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(im4, idx, p1, p2, p3, p4);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im4.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      if (Domain::IsInside(im4, i, j, this->_k, this->_l, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4);
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

    const int s1 =  im4.GetX() - (ei - bi);
    const int s2 = (im4.GetY() - (ej - bj)) * im4.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      if (Domain::IsInside(im4, i, j, k, this->_l, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4);
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
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  if (im4->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, vf);
  } else {
    QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  if (im4.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, vf);
  } else {
    QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  if (im4->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, vf);
  } else {
    QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
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
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  if (im4.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, vf);
  } else {
    QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
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
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, const GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, const GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// =============================================================================
// 3 const, 1 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 3 const, 1 non-const images
 */
template <class T1, class T2, class T3, class T4, class VoxelFunc>
struct QuaternaryForEachVoxelBody_3Const : public ForEachVoxelBody<VoxelFunc>
{
  const GenericImage<T1> &im1;
  const GenericImage<T2> &im2;
  const GenericImage<T3> &im3;
        GenericImage<T4> &im4;

  /// Constructor
  QuaternaryForEachVoxelBody_3Const(const GenericImage<T1> &im1,
                                    const GenericImage<T2> &im2,
                                    const GenericImage<T3> &im3,
                                          GenericImage<T4> &im4,
                                    VoxelFunc &vf)
  :
    ForEachVoxelBody<VoxelFunc>(vf, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4)
  {}

  /// Copy constructor
  QuaternaryForEachVoxelBody_3Const(const QuaternaryForEachVoxelBody_3Const &o)
  :
    ForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Split constructor
  QuaternaryForEachVoxelBody_3Const(QuaternaryForEachVoxelBody_3Const &o, split s)
  :
    ForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody_3Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody_3Const *>(this)->_VoxelFunc(im4, idx, p1, p2, p3, p4);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im4.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody_3Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4);
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

    const int s1 =  im4.GetX() - (ei - bi);
    const int s2 = (im4.GetY() - (ej - bj)) * im4.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody_3Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 3 const, 1 non-const images
 */
template <class T1, class T2, class T3, class T4,
          class VoxelFunc, class OutsideFunc = NaryVoxelFunction::NOP,
          class Domain = ForEachVoxelDomain::Foreground>
struct QuaternaryForEachVoxelIfBody_3Const : public ForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const GenericImage<T1> &im1;
  const GenericImage<T2> &im2;
  const GenericImage<T3> &im3;
        GenericImage<T4> &im4;

  /// Constructor
  QuaternaryForEachVoxelIfBody_3Const(const GenericImage<T1> &im1,
                                      const GenericImage<T2> &im2,
                                      const GenericImage<T3> &im3,
                                            GenericImage<T4> &im4,
                                      VoxelFunc &vf, OutsideFunc &of)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4)
  {}

  /// Copy constructor
  QuaternaryForEachVoxelIfBody_3Const(const QuaternaryForEachVoxelIfBody_3Const &o)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Split constructor
  QuaternaryForEachVoxelIfBody_3Const(QuaternaryForEachVoxelIfBody_3Const &o, split s)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4) {
      if (Domain::IsInside(im4, i, j, k, l, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody_3Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody_3Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      if (Domain::IsInside(im4, idx, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody_3Const *>(this)->_VoxelFunc  (im4, idx, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody_3Const *>(this)->_OutsideFunc(im4, idx, p1, p2, p3, p4);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im4.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      if (Domain::IsInside(im4, i, j, this->_k, this->_l, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody_3Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody_3Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4);
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

    const int s1 =  im4.GetX() - (ei - bi);
    const int s2 = (im4.GetY() - (ej - bj)) * im4.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      if (Domain::IsInside(im4, i, j, k, this->_l, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody_3Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody_3Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4);
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
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  if (im4->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, vf);
  } else {
    QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  if (im4.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, vf);
  } else {
    QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  if (im4->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, vf);
  } else {
    QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
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
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  if (im4.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, vf);
  } else {
    QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
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
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_3Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, const GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, const GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// =============================================================================
// 2 const, 2 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 2 const, 2 non-const images
 */
template <class T1, class T2, class T3, class T4, class VoxelFunc>
struct QuaternaryForEachVoxelBody_2Const : public ForEachVoxelBody<VoxelFunc>
{
  const GenericImage<T1> &im1;
  const GenericImage<T2> &im2;
        GenericImage<T3> &im3;
        GenericImage<T4> &im4;

  /// Constructor
  QuaternaryForEachVoxelBody_2Const(const GenericImage<T1> &im1,
                                    const GenericImage<T2> &im2,
                                          GenericImage<T3> &im3,
                                          GenericImage<T4> &im4,
                                    VoxelFunc &vf)
  :
    ForEachVoxelBody<VoxelFunc>(vf, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4)
  {}

  /// Copy constructor
  QuaternaryForEachVoxelBody_2Const(const QuaternaryForEachVoxelBody_2Const &o)
  :
    ForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Split constructor
  QuaternaryForEachVoxelBody_2Const(QuaternaryForEachVoxelBody_2Const &o, split s)
  :
    ForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody_2Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody_2Const *>(this)->_VoxelFunc(im4, idx, p1, p2, p3, p4);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im4.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody_2Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4);
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

    const int s1 =  im4.GetX() - (ei - bi);
    const int s2 = (im4.GetY() - (ej - bj)) * im4.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody_2Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 2 const, 2 non-const images
 */
template <class T1, class T2, class T3, class T4,
          class VoxelFunc, class OutsideFunc = NaryVoxelFunction::NOP,
          class Domain = ForEachVoxelDomain::Foreground>
struct QuaternaryForEachVoxelIfBody_2Const : public ForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const GenericImage<T1> &im1;
  const GenericImage<T2> &im2;
        GenericImage<T3> &im3;
        GenericImage<T4> &im4;

  /// Constructor
  QuaternaryForEachVoxelIfBody_2Const(const GenericImage<T1> &im1,
                                      const GenericImage<T2> &im2,
                                            GenericImage<T3> &im3,
                                            GenericImage<T4> &im4,
                                      VoxelFunc &vf, OutsideFunc &of)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4)
  {}

  /// Copy constructor
  QuaternaryForEachVoxelIfBody_2Const(const QuaternaryForEachVoxelIfBody_2Const &o)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Split constructor
  QuaternaryForEachVoxelIfBody_2Const(QuaternaryForEachVoxelIfBody_2Const &o, split s)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4) {
      if (Domain::IsInside(im4, i, j, k, l, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody_2Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody_2Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      if (Domain::IsInside(im4, idx, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody_2Const *>(this)->_VoxelFunc  (im4, idx, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody_2Const *>(this)->_OutsideFunc(im4, idx, p1, p2, p3, p4);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im4.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      if (Domain::IsInside(im4, i, j, this->_k, this->_l, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody_2Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody_2Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4);
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

    const int s1 =  im4.GetX() - (ei - bi);
    const int s2 = (im4.GetY() - (ej - bj)) * im4.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      if (Domain::IsInside(im4, i, j, k, this->_l, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody_2Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody_2Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4);
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
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  if (im4->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, vf);
  } else {
    QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  if (im4.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, vf);
  } else {
    QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  if (im4->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, vf);
  } else {
    QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
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
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  if (im4.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, vf);
  } else {
    QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
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
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_2Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, const GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, const GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// =============================================================================
// 1 const, 3 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 1 const, 3 non-const images
 */
template <class T1, class T2, class T3, class T4, class VoxelFunc>
struct QuaternaryForEachVoxelBody_1Const : public ForEachVoxelBody<VoxelFunc>
{
  const GenericImage<T1> &im1;
        GenericImage<T2> &im2;
        GenericImage<T3> &im3;
        GenericImage<T4> &im4;

  /// Constructor
  QuaternaryForEachVoxelBody_1Const(const GenericImage<T1> &im1,
                                          GenericImage<T2> &im2,
                                          GenericImage<T3> &im3,
                                          GenericImage<T4> &im4,
                                    VoxelFunc &vf)
  :
    ForEachVoxelBody<VoxelFunc>(vf, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4)
  {}

  /// Copy constructor
  QuaternaryForEachVoxelBody_1Const(const QuaternaryForEachVoxelBody_1Const &o)
  :
    ForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Split constructor
  QuaternaryForEachVoxelBody_1Const(QuaternaryForEachVoxelBody_1Const &o, split s)
  :
    ForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(im4, idx, p1, p2, p3, p4);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im4.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4);
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

    const int s1 =  im4.GetX() - (ei - bi);
    const int s2 = (im4.GetY() - (ej - bj)) * im4.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 1 const, 3 non-const images
 */
template <class T1, class T2, class T3, class T4,
          class VoxelFunc, class OutsideFunc = NaryVoxelFunction::NOP,
          class Domain = ForEachVoxelDomain::Foreground>
struct QuaternaryForEachVoxelIfBody_1Const : public ForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const GenericImage<T1> &im1;
        GenericImage<T2> &im2;
        GenericImage<T3> &im3;
        GenericImage<T4> &im4;

  /// Constructor
  QuaternaryForEachVoxelIfBody_1Const(const GenericImage<T1> &im1,
                                            GenericImage<T2> &im2,
                                            GenericImage<T3> &im3,
                                            GenericImage<T4> &im4,
                                      VoxelFunc &vf, OutsideFunc &of)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4)
  {}

  /// Copy constructor
  QuaternaryForEachVoxelIfBody_1Const(const QuaternaryForEachVoxelIfBody_1Const &o)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Split constructor
  QuaternaryForEachVoxelIfBody_1Const(QuaternaryForEachVoxelIfBody_1Const &o, split s)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4) {
      if (Domain::IsInside(im4, i, j, k, l, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      if (Domain::IsInside(im4, idx, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (im4, idx, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(im4, idx, p1, p2, p3, p4);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im4.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      if (Domain::IsInside(im4, i, j, this->_k, this->_l, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4);
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

    const int s1 =  im4.GetX() - (ei - bi);
    const int s2 = (im4.GetY() - (ej - bj)) * im4.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      if (Domain::IsInside(im4, i, j, k, this->_l, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4);
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
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  if (im4->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, vf);
  } else {
    QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  if (im4.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, vf);
  } else {
    QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  if (im4->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, vf);
  } else {
    QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
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
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  if (im4.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, vf);
  } else {
    QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
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
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody_1Const<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// =============================================================================
// 4 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 4 non-const images
 */
template <class T1, class T2, class T3, class T4, class VoxelFunc>
struct QuaternaryForEachVoxelBody : public ForEachVoxelBody<VoxelFunc>
{
  GenericImage<T1> &im1;
  GenericImage<T2> &im2;
  GenericImage<T3> &im3;
  GenericImage<T4> &im4;

  /// Constructor
  QuaternaryForEachVoxelBody(GenericImage<T1> &im1,
                             GenericImage<T2> &im2,
                             GenericImage<T3> &im3,
                             GenericImage<T4> &im4,
                             VoxelFunc &vf)
  :
    ForEachVoxelBody<VoxelFunc>(vf, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4)
  {}

  /// Copy constructor
  QuaternaryForEachVoxelBody(const QuaternaryForEachVoxelBody &o)
  :
    ForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Split constructor
  QuaternaryForEachVoxelBody(QuaternaryForEachVoxelBody &o, split s)
  :
    ForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody *>(this)->_VoxelFunc(im4, idx, p1, p2, p3, p4);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im4.GetX() - (ei - bi);

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4);
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

    const int s1 =  im4.GetX() - (ei - bi);
    const int s2 = (im4.GetY() - (ej - bj)) * im4.GetX();

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<QuaternaryForEachVoxelBody *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 4 non-const images
 */
template <class T1, class T2, class T3, class T4,
          class VoxelFunc, class OutsideFunc = NaryVoxelFunction::NOP,
          class Domain = ForEachVoxelDomain::Foreground>
struct QuaternaryForEachVoxelIfBody : public ForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  GenericImage<T1> &im1;
  GenericImage<T2> &im2;
  GenericImage<T3> &im3;
  GenericImage<T4> &im4;

  /// Constructor
  QuaternaryForEachVoxelIfBody(GenericImage<T1> &im1,
                               GenericImage<T2> &im2,
                               GenericImage<T3> &im3,
                               GenericImage<T4> &im4,
                               VoxelFunc &vf, OutsideFunc &of)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of, im1.Attributes()), im1(im1), im2(im2), im3(im3), im4(im4)
  {}

  /// Copy constructor
  QuaternaryForEachVoxelIfBody(const QuaternaryForEachVoxelIfBody &o)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Split constructor
  QuaternaryForEachVoxelIfBody(QuaternaryForEachVoxelIfBody &o, split s)
  :
    ForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4)
  {}

  /// Process entire image
  void operator ()(const ImageAttributes &attr) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4) {
      if (Domain::IsInside(im4, i, j, k, l, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4);
    }
  }

  /// Process image region using linear index
  void operator ()(const blocked_range<int> &re) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      if (Domain::IsInside(im4, idx, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody *>(this)->_VoxelFunc  (im4, idx, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody *>(this)->_OutsideFunc(im4, idx, p1, p2, p3, p4);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im4.GetX() - (ei - bi);

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      if (Domain::IsInside(im4, i, j, this->_k, this->_l, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4);
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

    const int s1 =  im4.GetX() - (ei - bi);
    const int s2 = (im4.GetY() - (ej - bj)) * im4.GetX();

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1) {
      if (Domain::IsInside(im4, i, j, k, this->_l, p4)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<QuaternaryForEachVoxelIfBody *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4);
      } else const_cast<QuaternaryForEachVoxelIfBody *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4);
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
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  if (im4->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, vf);
  } else {
    QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  if (im4.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, vf);
  } else {
    QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  if (im4->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, vf);
  } else {
    QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
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
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(*im1, *im2, *im3, *im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  if (im4.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, vf);
  } else {
    QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
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
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  QuaternaryForEachVoxelBody<T1, T2, T3, T4, VoxelFunc> body(im1, im2, im3, im4, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  blocked_range<int> re(0, im4->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
    blocked_range<int> re(0, im4->GetNumberOfVoxels() / im4->GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, GenericImage<T1> *im1, GenericImage<T2> *im2, GenericImage<T3> *im3, GenericImage<T4> *im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  blocked_range<int> re(0, im4.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  if (im4.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
  } else {
    QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
    blocked_range<int> re(0, im4.GetNumberOfVoxels() / im4.GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const ImageAttributes &attr, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf, OutsideFunc &of)
{
  QuaternaryForEachVoxelIfBody<T1, T2, T3, T4, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4, VoxelFunc &vf)
{
  NaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, GenericImage<T1> &im1, GenericImage<T2> &im2, GenericImage<T3> &im3, GenericImage<T4> &im4)
{
  if (VoxelFunc::IsReduction()) _foreachquaternaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, vf);
}


} // namespace mirtk

#endif
