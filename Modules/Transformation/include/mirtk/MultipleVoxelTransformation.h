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
 */

#ifndef MIRTK_MultipleVoxelTransformation_H
#define MIRTK_MultipleVoxelTransformation_H

#include "mirtk/Float.h"
#include "mirtk/GenericImage.h"
#include "mirtk/VoxelFunction.h"
#include "mirtk/InterpolateImageFunction.h"

#include <cstdlib>
#include <limits>
#include <iostream>


namespace mirtk {


/**
 * These voxel functions are used by MultipleImageTransformation
 * to efficiently transform one or more input images at once
 */
namespace MultipleVoxelTransformation {


////////////////////////////////////////////////////////////////////////////////
// Base class for voxel transformation functions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
/// Base class for voxel transformation functions
class Base : public VoxelFunction
{
protected:

  /// Default constructor
  Base(int numvox = 0)
  :
    _NumberOfVoxels     (numvox),
    _TwiceNumberOfVoxels(numvox<<1)
  {}

  /// Copy constructor
  Base(const Base &other)
  :
    _NumberOfVoxels     (other._NumberOfVoxels),
    _TwiceNumberOfVoxels(other._TwiceNumberOfVoxels)
  {}

  int _NumberOfVoxels;      ///< Number of voxels per vector component
  int _TwiceNumberOfVoxels; ///< Offset of z component

public:

  /// Get reference to x component of vector field
  template <class T>
  inline       T &x_(      T *p) const { return p[0]; }

  /// Get const reference to x component of vector field
  template <class T>
  inline const T &x_(const T *p) const { return p[0]; }

  /// Get reference to y component of vector field
  template <class T>
  inline       T &y_(      T *p) const { return p[_NumberOfVoxels]; }

  /// Get const reference to y component of vector field
  template <class T>
  inline const T &y_(const T *p) const { return p[_NumberOfVoxels]; }

  /// Get reference to z component of vector field
  template <class T>
  inline       T &z_(      T *p) const { return p[_TwiceNumberOfVoxels]; }

  /// Get const reference to z component of vector field
  template <class T>
  inline const T &z_(const T *p) const { return p[_TwiceNumberOfVoxels]; }
};

// -----------------------------------------------------------------------------
/// Input/output data members of voxel transformation functions
struct TransformData
{
  int                        _NumberOfImages;        ///< Number of input images
  int                        _MaxNumberOfComponents; ///< Maximum number of image components/channels (_t)
  const BaseImage          **_Inputs;                ///< Input images
  const BaseImage           *_Input;                 ///< Single/reference input
  InterpolateImageFunction **_Interpolators;         ///< Input interpolators
  InterpolateImageFunction  *_Interpolator;          ///< Single/reference interpolator
  double                    *_PaddingValues;         ///< Outside padding values
  double                     _PaddingValue;          ///< Default padding value
  const Transformation      *_Transformation1;       ///< First transformation
  const Transformation      *_Transformation2;       ///< Second transformation (optional)
  bool                       _Invert;                ///< Whether to invert the transformations
  double                     _ScaleFactor;           ///< Output intensity scale factor
  double                     _Offset;                ///< Output intensity offset
  Image                    **_Outputs;               ///< Output images
  Image                     *_Output;                ///< Single/reference output image
  BinaryImage               *_Mask;                  ///< Output foreground mask
  int                        _Begin;                 ///< Index of first image to process
  int                        _End;                   ///< Index of last  image to process + 1

  /// Default constructor
  TransformData()
  :
    _NumberOfImages       (0),
    _MaxNumberOfComponents(0),
    _Inputs               (NULL),
    _Input                (NULL),
    _Interpolators        (NULL),
    _Interpolator         (NULL),
    _PaddingValues        (NULL),
    _PaddingValue         (.0),
    _Transformation1      (NULL),
    _Transformation2      (NULL),
    _Invert               (false),
    _ScaleFactor          (1.0),
    _Offset               (0.0),
    _Outputs              (NULL),
    _Output               (NULL),
    _Mask                 (NULL),
    _Begin                (-1),
    _End                  (-1)
  {}
};

// -----------------------------------------------------------------------------
/**
 * Base class for voxel transformation functions with various transformation methods
 *
 * Using an explicit interpolate image function type helps the compiler to deduce
 * which virtual methods are being called and thus inline the code for better
 * performance. If the actual type of the interpolator is not known when
 * instantiating this class, use the default InterpolateImageFunction base
 * interpolate image function type.
 */
template <class InterpolateImageFunction = InterpolateImageFunction,
          class InputDomain              = InterpolationDomain::Foreground>
struct BaseTransform : public TransformData, public Base
{
  // ---------------------------------------------------------------------------
  /// Get (casted) pointer to interpolator of type InterpolateImageFunction
  InterpolateImageFunction *GetInterpolator(int n)
  {
    return reinterpret_cast<InterpolateImageFunction *>(_Interpolators[n]);
  }

  // ---------------------------------------------------------------------------
  /// Get (casted) pointer to reference interpolator of type InterpolateImageFunction
  InterpolateImageFunction *GetInterpolator()
  {
    return reinterpret_cast<InterpolateImageFunction *>(_Interpolator);
  }

  // ===========================================================================
  // Construction/Destruction
  // ===========================================================================

protected:

  // -------------------------------------------------------------------------
  /// Constructor
  BaseTransform(const TransformData &data)
  :
    TransformData(data), Base()
  {
    // Allow user to specify only the main input/output if only one image to process
    if (!_Inputs || !_Interpolators || !_Outputs) {
      if (_NumberOfImages > 1) {
        cerr << "MultipleVoxelTransformation::BaseTransform: Invalid TransformData" << endl;
        exit(1);
      }
      if (!_Inputs) {
        if (!_Input) {
          cerr << "MultipleVoxelTransformation::BaseTransform: Missing input image" << endl;
          exit(1);
        }
        _Inputs = &_Input;
      }
      if (!_Interpolators) {
        if (!_Interpolator) {
          cerr << "MultipleVoxelTransformation::BaseTransform: Missing input interpolator" << endl;
          exit(1);
        }
        _Interpolators = &_Interpolator;
      }
      if (!_Outputs) {
        if (!_Output) {
          cerr << "MultipleVoxelTransformation::BaseTransform: Missing output image" << endl;
          exit(1);
        }
        _Outputs = &_Output;
      }
      _NumberOfImages = 1;
    }
    // Set main input/output if not set yet
    if (!_Input)        _Input        = _Inputs       [0];
    if (!_Interpolator) _Interpolator = _Interpolators[0];
    if (!_Output)       _Output       = _Outputs      [0];
    // Set offsets for vector component access
    _NumberOfVoxels      = _Output->GetX() * _Output->GetY() * _Output->GetZ();
    _TwiceNumberOfVoxels = 2 * _NumberOfVoxels;
    // Range of output images to process
    if (_Begin < 0) { _Begin = 0;  _End = -1; }
    if (_End   < 0) /*_Begin = x*/ _End = _NumberOfImages;
    if (_End   < _Begin) _End = _Begin;
    // Determine maximum number of components
    if (_MaxNumberOfComponents < 1) {
      _MaxNumberOfComponents = _Inputs[0]->GetT();
      for (int n = 1; n < _NumberOfImages; ++n) {
        _MaxNumberOfComponents = max(_MaxNumberOfComponents, _Inputs[n]->GetT());
      }
    }
    // Allocate memory for interpolated values
    _v = new double[_MaxNumberOfComponents];
  }

  // -------------------------------------------------------------------------
  /// Copy constructor
  BaseTransform(const BaseTransform &other)
  :
    TransformData(other), Base(other)
  {
    _v = new double[_MaxNumberOfComponents];
  }

public:

  // ---------------------------------------------------------------------------
  /// Destructor
  virtual ~BaseTransform()
  {
    delete[] _v;
  }

  // ===========================================================================
  // 1. Output voxel indices to world coordinates
  // ===========================================================================

  // ---------------------------------------------------------------------------
  /// Convert output voxel indices to world coordinates
  inline void OutputToWorld(int i, int j, int k)
  {
    _x = static_cast<double>(i);
    _y = static_cast<double>(j);
    _z = static_cast<double>(k);
    _Output->ImageToWorld(_x, _y, _z);
  }

  // ---------------------------------------------------------------------------
  /// Convert output voxel indices to world coordinates given pointer to lookup table
  inline void OutputToWorld(const double *i2w)
  {
    _x = x_(i2w);
    _y = y_(i2w);
    _z = z_(i2w);
  }

  // ===========================================================================
  // 2. World coordinate transformation
  // ===========================================================================

  // ---------------------------------------------------------------------------
  /// Transform indices using single transformation
  ///
  /// Ignores the second transformation assuming it to be unused.
  inline void ApplyTransformation()
  {
    if (_Invert) _Transformation1->Inverse  (_x, _y, _z);
    else         _Transformation1->Transform(_x, _y, _z);
  }

  // ---------------------------------------------------------------------------
  /// Transform indices using single dense displacement field
  ///
  /// Requires dense displacement field computed from _Transformation1 and
  /// ignores the second transformation, assuming it to be unused.
  inline void ApplyDisplacement(const double *disp1)
  {
    _x += x_(disp1);
    _y += y_(disp1);
    _z += z_(disp1);
  }

  // ---------------------------------------------------------------------------
  /// Transform using fluid composition of transformations
  inline void ApplyTransformations()
  {
    if (_Invert) {
      _Transformation1->Inverse  (_x, _y, _z);
      _Transformation2->Inverse  (_x, _y, _z);
    } else {
      _Transformation1->Transform(_x, _y, _z);
      _Transformation2->Transform(_x, _y, _z);
    }
  }

  // ---------------------------------------------------------------------------
  /// Transform indices using fluid composition of dense displacement field and transformation
  ///
  /// Requires dense displacement field computed from _Transformation1.
  inline void ApplyDisplacementAndTransformation(const double *disp1)
  {
    _x += x_(disp1);
    _y += y_(disp1);
    _z += z_(disp1);
    if (_Invert) _Transformation2->Inverse  (_x, _y, _z);
    else         _Transformation2->Transform(_x, _y, _z);
  }

  // ---------------------------------------------------------------------------
  /// Transform indices using additive composition of dense displacement fields
  ///
  /// Requires dense displacement fields computed from _Transformation1 and _Transformation2.
  inline void ApplyDisplacements(const double *disp1, const double *disp2)
  {
    _x += x_(disp1) + x_(disp2);
    _y += y_(disp1) + y_(disp2);
    _z += z_(disp1) + z_(disp2);
  }

  // ===========================================================================
  // 3. World coordinates to input voxel indices
  // ===========================================================================

  // ---------------------------------------------------------------------------
  /// Convert world coordinates to input voxel indices
  inline void WorldToInput()
  {
    _Input->WorldToImage(_x, _y, _z);
  }

  // ===========================================================================
  // 4. i) Either store voxel transformation in dense vector field...
  //
  // Note: Faster if output images have varying scalar type as then each output
  //       can be updated separately using the pre-computed voxel index
  //       index transformation map and an update voxel function which is
  //       specialized for each output's particular scalar type.
  // ===========================================================================

  // ---------------------------------------------------------------------------
  /// Update map from output to input voxel indices with inside check
  inline void PutVoxelTransformation(int i, int j, int k, double *o2i)
  {
    if (InputDomain::IsInside(GetInterpolator(), _x, _y, _z)) {
      x_(o2i) = _x;
      y_(o2i) = _y;
      z_(o2i) = _z;
      if (_Mask) _Mask->Put(i, j, k, static_cast<BinaryPixel>(true));
    } else {
      x_(o2i) = numeric_limits<double>::quiet_NaN();
      if (_Mask) _Mask->Put(i, j, k, static_cast<BinaryPixel>(false));
    }
  }

  // ===========================================================================
  // 4. ii) ...or set directly all output voxels at once to interpolated values
  //
  // Note: Faster and more memory efficient when all output images have the
  //       same scalar type as no intermediate storage and lookup of the
  //       transformed voxel coordinates is required. This is in particular the
  //       case when only a single input image is transformed.
  // ===========================================================================

  // ---------------------------------------------------------------------------
  /// Set outputs to dedicated outside value
  ///
  /// Efficient and applicable when all output images have a known unique scalar type.
  template <class OutputVoxelType>
  inline void PutOutsideValue(int i, int j, int k)
  {
    // Set output voxel to outside value
    OutputVoxelType *p;
    for (int n = _Begin; n != _End; ++n) {
      p = reinterpret_cast<OutputVoxelType *>(_Outputs[n]->GetScalarPointer(i, j, k));
      for (int l = 0; l < _Outputs[n]->GetT(); ++l, p += _NumberOfVoxels) {
        *p = static_cast<OutputVoxelType>(_PaddingValues ? _PaddingValues[n] : _PaddingValue);
      }
    }
    // Update foreground mask
    if (_Mask) _Mask->Put(i, j, k, static_cast<BinaryPixel>(false));
  }

  // ---------------------------------------------------------------------------
  /// Set outputs to dedicated outside value
  ///
  /// Less efficient than PutOutsideValue but can also be used when output
  /// images have not a known unique scalar type.
  inline void PutOutsideValueAsDouble(int i, int j, int k)
  {
    // Set output voxel to outside value
    for (int n = _Begin; n != _End; ++n) {
      for (int l = 0; l < _Outputs[n]->GetT(); ++l) {
        _Outputs[n]->PutAsDouble(i, j, k, l, _PaddingValues ? _PaddingValues[n] : _PaddingValue);
      }
    }
    // Update foreground mask
    if (_Mask) _Mask->Put(i, j, k, static_cast<BinaryPixel>(false));
  }

  // ---------------------------------------------------------------------------
  /// Set outputs to interpolated input values
  ///
  /// Most efficient implementation. Only applicable when all output images have
  /// double as scalar type and thus also no intensity rescaling is required.
  inline void InterpolatePut(int i, int j, int k, double *)
  {
    if (InputDomain::IsInside(GetInterpolator(), _x, _y, _z)) {
      // Set to interpolated value
      for (int n = _Begin; n != _End; ++n) {
        GetInterpolator(n)->EvaluateInside(
            reinterpret_cast<double *>(_Outputs[n]->GetScalarPointer(i, j, k)), _x, _y, _z, _NumberOfVoxels);
      }
      // Update foreground mask
      if (_Mask) _Mask->Put(i, j, k, static_cast<BinaryPixel>(true));
    } else {
      // Set to outside value
      PutOutsideValue<double>(i, j, k);
    }
  }

  // ---------------------------------------------------------------------------
  /// Set outputs to interpolated input values
  ///
  /// Slightly less efficient implementation. Applicable when all output images
  /// have a known unique scalar type and no intensity rescaling is used.
  template <class OutputVoxelType>
  inline void InterpolatePut(int i, int j, int k, OutputVoxelType *)
  {
    if (InputDomain::IsInside(GetInterpolator(), _x, _y, _z)) {
      // Set to interpolated value
      OutputVoxelType *p;
      for (int n = _Begin; n != _End; ++n) {
        GetInterpolator(n)->EvaluateInside(_v, _x, _y, _z);
        p = reinterpret_cast<OutputVoxelType *>(_Outputs[n]->GetScalarPointer(i, j, k));
        for (int l = 0; l < _Outputs[n]->GetT(); ++l, p += _NumberOfVoxels) {
          *p = static_cast<OutputVoxelType>(_v[l]);
        }
      }
      // Update foreground mask
      if (_Mask) _Mask->Put(i, j, k, static_cast<BinaryPixel>(true));
    } else {
      // Set to outside value
      PutOutsideValue<OutputVoxelType>(i, j, k);
    }
  }

  // ---------------------------------------------------------------------------
  /// Set outputs to interpolated input values
  ///
  /// Least efficient method without rescaling. Does not imply any particular
  /// output scalar type or interpolator type. In this case it is better to
  /// first store the voxel index transformation in the _Output2Input map and
  /// then update each output image separately.
  inline void InterpolatePutAsDouble(int i, int j, int k)
  {
    if (InputDomain::IsInside(GetInterpolator(), _x, _y, _z)) {
      // Set to interpolated value
      for (int n = _Begin; n != _End; ++n) {
        GetInterpolator(n)->EvaluateInside(_v, _x, _y, _z);
        for (int l = 0; l < _Outputs[n]->GetT(); ++l) {
          _Outputs[n]->PutAsDouble(i, j, k, l, _v[k]);
        }
      }
      // Update foreground mask
      if (_Mask) _Mask->Put(i, j, k, static_cast<BinaryPixel>(true));
    } else {
      // Set to outside value
      PutOutsideValueAsDouble(i, j, k);
    }
  }

  // ---------------------------------------------------------------------------
  /// Set outputs to interpolated and rescaled input values
  ///
  /// Efficient and applicable when all output images have a known unique scalar type.
  template <class OutputVoxelType>
  inline void InterpolateRescalePut(int i, int j, int k, OutputVoxelType *)
  {
    if (InputDomain::IsInside(GetInterpolator(), _x, _y, _z)) {
      // Set to interpolated and rescaled value
      OutputVoxelType *p;
      for (int n = _Begin; n != _End; ++n) {
        GetInterpolator(n)->EvaluateInside(_v, _x, _y, _z);
        p = reinterpret_cast<OutputVoxelType *>(_Outputs[n]->GetScalarPointer(i, j, k));
        for (int l = 0; l < _Outputs[n]->GetT(); ++l, p += _NumberOfVoxels) {
          *p = static_cast<OutputVoxelType>(_ScaleFactor * _v[l] + _Offset);
        }
      }
      // Update foreground mask
      if (_Mask) _Mask->Put(i, j, k, static_cast<BinaryPixel>(true));
    } else {
      // Set to outside value
      PutOutsideValue<OutputVoxelType>(i, j, k);
    }
  }

  // ---------------------------------------------------------------------------
  /// Set outputs to interpolated and rescaled input values
  ///
  /// Least efficient but most general interpolation method. Does not imply any
  /// particular output scalar type or interpolator type. In this case it is
  /// better, however, to first store the voxel index transformation in the
  /// _Output2Input map and then update each output image separately if the
  /// additional required memory is not an issue.
  inline void InterpolateRescalePutAsDouble(int i, int j, int k)
  {
    if (InputDomain::IsInside(GetInterpolator(), _x, _y, _z)) {
      // Set to interpolated and rescaled value
      for (int n = _Begin; n != _End; ++n) {
        GetInterpolator(n)->EvaluateInside(_v, _x, _y, _z);
        for (int l = 0; l < _Outputs[n]->GetT(); ++l) {
          _Outputs[n]->PutAsDouble(i, j, k, l, _ScaleFactor * _v[l] + _Offset);
        }
      }
      // Update foreground mask
      if (_Mask) _Mask->Put(i, j, k, static_cast<BinaryPixel>(true));
    } else {
      // Set to outside value
      PutOutsideValueAsDouble(i, j, k);
    }
  }

private:

  double  _x, _y, _z; ///< (Intermediate) Transformed world/image coordinates
  double *_v;         ///< Pre-allocated memory for interpolated values

};

////////////////////////////////////////////////////////////////////////////////
// Voxel transformation functions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct MapIndices : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  MapIndices(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  void operator ()(int i, int j, int k, int, double *o2i)
  {
    this->OutputToWorld(i, j, k);
    this->WorldToInput();
    this->PutVoxelTransformation(i, j, k, o2i);
  }

  void operator ()(int i, int j, int k, int, const double *o2w, double *o2i)
  {
    this->OutputToWorld(o2w);
    this->WorldToInput();
    this->PutVoxelTransformation(i, j, k, o2i);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct Transform : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  Transform(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  void operator ()(int i, int j, int k, int, double *o2i)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyTransformation();
    this->WorldToInput();
    this->PutVoxelTransformation(i, j, k, o2i);
  }

  void operator ()(int i, int j, int k, int, const double *o2w, double *o2i)
  {
    this->OutputToWorld(o2w);
    this->ApplyTransformation();
    this->WorldToInput();
    this->PutVoxelTransformation(i, j, k, o2i);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct Displace : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  Displace(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  void operator ()(int i, int j, int k, int, const double *disp1, double *o2i)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyDisplacement(disp1);
    this->WorldToInput();
    this->PutVoxelTransformation(i, j, k, o2i);
  }

  void operator ()(int i, int j, int k, int, const double *o2w, const double *disp1, double *o2i)
  {
    this->OutputToWorld(o2w);
    this->ApplyDisplacement(disp1);
    this->WorldToInput();
    this->PutVoxelTransformation(i, j, k, o2i);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct FluidTransform : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  FluidTransform(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  void operator ()(int i, int j, int k, int, double *o2i)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyTransformations();
    this->WorldToInput();
    this->PutVoxelTransformation(i, j, k, o2i);
  }

  void operator ()(int i, int j, int k, int, const double *o2w, double *o2i)
  {
    this->OutputToWorld(o2w);
    this->ApplyTransformations();
    this->WorldToInput();
    this->PutVoxelTransformation(i, j, k, o2i);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct DisplaceTransform : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  DisplaceTransform(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  void operator ()(int i, int j, int k, int, const double *disp1, double *o2i)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyDisplacementAndTransformation(disp1);
    this->WorldToInput();
    this->PutVoxelTransformation(i, j, k, o2i);
  }

  void operator ()(int i, int j, int k, int, const double *o2w, const double *disp1, double *o2i)
  {
    this->OutputToWorld(o2w);
    this->ApplyDisplacementAndTransformation(disp1);
    this->WorldToInput();
    this->PutVoxelTransformation(i, j, k, o2i);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct AddDisplacements : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  AddDisplacements(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  void operator ()(int i, int j, int k, int, const double *disp1, const double *disp2, double *o2i)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyDisplacements(disp1, disp2);
    this->WorldToInput();
    this->PutVoxelTransformation(i, j, k, o2i);
  }

  void operator ()(int i, int j, int k, int, const double *o2w, const double *disp1, const double *disp2, double *o2i)
  {
    this->OutputToWorld(o2w);
    this->ApplyDisplacements(disp1, disp2);
    this->WorldToInput();
    this->PutVoxelTransformation(i, j, k, o2i);
  }
};

////////////////////////////////////////////////////////////////////////////////
// Voxel interpolation functions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
/// Set output image of known scalar type to interpolated value using pre-computed
/// map from output voxel indices to input indices
template <class InterpolateImageFunction = InterpolateImageFunction>
struct Interpolate : public Base
{
  InterpolateImageFunction *_Interpolator;
  int                       _NumberOfComponents;
  double                    _PaddingValue;

  /// Default constructor
  Interpolate() : Base(), _v(NULL) {}

  /// Copy constructor
  Interpolate(const Interpolate &) : _v(NULL) {}

  /// Destructor
  virtual ~Interpolate()
  {
    delete[] _v;
  }

  /// Allocate memory for interpolated values the first time
  inline void InitializeMemory()
  {
    if (!_v) _v = new double[_NumberOfComponents];
  }

  template <class OutputVoxelType>
  void operator ()(int, int, int, int, const double *o2i, OutputVoxelType *out)
  {
    if (IsNaN(*o2i)) {
      for (int l = 0; l < _NumberOfComponents; ++l, out += _NumberOfVoxels) {
        *out = static_cast<OutputVoxelType>(_PaddingValue);
      }
    } else {
      if (!_v) _v = new double[_NumberOfComponents];
      _Interpolator->EvaluateInside(_v, x_(o2i), y_(o2i), z_(o2i));
      for (int l = 0; l < _NumberOfComponents; ++l, out += _NumberOfVoxels) {
        *out = static_cast<OutputVoxelType>(_v[l]);
      }
    }
  }

  void operator ()(int, int, int, int, const double *o2i, double *out)
  {
    if (IsNaN(*o2i)) {
      for (int l = 0; l < _NumberOfComponents; ++l, out += _NumberOfVoxels) {
        *out = _PaddingValue;
      }
    } else {
      _Interpolator->EvaluateInside(out, x_(o2i), y_(o2i), z_(o2i), _NumberOfVoxels);
    }
  }

protected:
  double *_v; ///< Pre-allocated memory for interpolated values
};

// -----------------------------------------------------------------------------
/// Set output image of known scalar type to interpolated and rescaled value
/// using pre-computed map from output voxel indices to input indices
template <class InterpolateImageFunction = InterpolateImageFunction>
struct InterpolateRescale : public Interpolate<InterpolateImageFunction>
{
  double _ScaleFactor;
  double _Offset;

  template <class OutputVoxelType>
  void operator ()(int, int, int, int, const double *o2i, OutputVoxelType *out)
  {
    if (IsNaN(*o2i)) {
      for (int l = 0; l < this->_NumberOfComponents; ++l, out += this->_NumberOfVoxels) {
        *out = static_cast<OutputVoxelType>(this->_PaddingValue);
      }
    } else {
      this->InitializeMemory();
      this->_Interpolator->EvaluateInside(this->_v, this->x_(o2i), this->y_(o2i), this->z_(o2i));
      for (int l = 0; l < this->_NumberOfComponents; ++l, out += this->_NumberOfVoxels) {
        *out = static_cast<OutputVoxelType>(_ScaleFactor * this->_v[l] + _Offset);
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Set output image of unknown scalar type to interpolated and rescaled value
/// using pre-computed map from output voxel indices to input indices
template <class InterpolateImageFunction = InterpolateImageFunction>
struct InterpolateRescaleAsDouble : public InterpolateRescale<InterpolateImageFunction>
{
  Image *_Output;

  void operator ()(int i, int j, int k, int, const double *o2i)
  {
    if (IsNaN(*o2i)) {
      for (int l = 0; l < this->_NumberOfComponents; ++l) {
        this->_Output->PutAsDouble(i, j, k, l, this->_PaddingValue);
      }
    } else {
      this->InitializeMemory();
      this->_Interpolator->EvaluateInside(this->_v, this->x_(o2i), this->y_(o2i), this->z_(o2i));
      for (int l = 0; l < this->_NumberOfComponents; ++l) {
        _Output->PutAsDouble(i, j, k, l, this->_ScaleFactor * this->_v[l] + this->_Offset);
      }
    }
  }
};

////////////////////////////////////////////////////////////////////////////////
// Voxel transformation and interpolation functions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
// Apply world coordinates
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct MapIndicesInterpolatePut : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  MapIndicesInterpolatePut(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, OutputVoxelType *out)
  {
    this->OutputToWorld(i, j, k);
    this->WorldToInput();
    this->template InterpolatePut(i, j, k, out);
  }

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, const double *o2w, OutputVoxelType *out)
  {
    this->OutputToWorld(o2w);
    this->WorldToInput();
    this->template InterpolatePut(i, j, k, out);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct MapIndicesInterpolatePutAsDouble : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  MapIndicesInterpolatePutAsDouble(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  void operator ()(int i, int j, int k, int, void *)
  {
    this->OutputToWorld(i, j, k);
    this->WorldToInput();
    this->InterpolatePutAsDouble(i, j, k);
  }

  void operator ()(int i, int j, int k, int, const double *o2w, void *)
  {
    this->OutputToWorld(o2w);
    this->WorldToInput();
    this->InterpolatePutAsDouble(i, j, k);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct MapIndicesInterpolateRescalePut : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  MapIndicesInterpolateRescalePut(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, OutputVoxelType *out)
  {
    this->OutputToWorld(i, j, k);
    this->WorldToInput();
    this->template InterpolateRescalePut(i, j, k, out);
  }

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, const double *o2w, OutputVoxelType *out)
  {
    this->OutputToWorld(o2w);
    this->WorldToInput();
    this->template InterpolateRescalePut(i, j, k, out);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct MapIndicesInterpolateRescalePutAsDouble : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  MapIndicesInterpolateRescalePutAsDouble(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  void operator ()(int i, int j, int k, int, void *)
  {
    this->OutputToWorld(i, j, k);
    this->WorldToInput();
    this->InterpolateRescalePutAsDouble(i, j, k);
  }

  void operator ()(int i, int j, int k, int, const double *o2w, void *)
  {
    this->OutputToWorld(o2w);
    this->WorldToInput();
    this->InterpolateRescalePutAsDouble(i, j, k);
  }
};

// -----------------------------------------------------------------------------
// Apply single transformation
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct TransformInterpolatePut : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  TransformInterpolatePut(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, OutputVoxelType *out)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyTransformation();
    this->WorldToInput();
    this->template InterpolatePut(i, j, k, out);
  }

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, const double *o2w, OutputVoxelType *out)
  {
    this->OutputToWorld(o2w);
    this->ApplyTransformation();
    this->WorldToInput();
    this->template InterpolatePut(i, j, k, out);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct TransformInterpolatePutAsDouble : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  TransformInterpolatePutAsDouble(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  void operator ()(int i, int j, int k, int, void *)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyTransformation();
    this->WorldToInput();
    this->InterpolatePutAsDouble(i, j, k);
  }

  void operator ()(int i, int j, int k, int, const double *o2w, void *)
  {
    this->OutputToWorld(o2w);
    this->ApplyTransformation();
    this->WorldToInput();
    this->InterpolatePutAsDouble(i, j, k);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct TransformInterpolateRescalePut : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  TransformInterpolateRescalePut(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, OutputVoxelType *out)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyTransformation();
    this->WorldToInput();
    this->template InterpolateRescalePut(i, j, k, out);
  }

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, const double *o2w, OutputVoxelType *out)
  {
    this->OutputToWorld(o2w);
    this->ApplyTransformation();
    this->WorldToInput();
    this->template InterpolateRescalePut(i, j, k, out);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct TransformInterpolateRescalePutAsDouble : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  TransformInterpolateRescalePutAsDouble(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  void operator ()(int i, int j, int k, int, void *)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyTransformation();
    this->WorldToInput();
    this->InterpolateRescalePutAsDouble(i, j, k);
  }

  void operator ()(int i, int j, int k, int, const double *o2w, void *)
  {
    this->OutputToWorld(o2w);
    this->ApplyTransformation();
    this->WorldToInput();
    this->InterpolateRescalePutAsDouble(i, j, k);
  }
};

// -----------------------------------------------------------------------------
// Apply single displacement
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct DisplaceInterpolatePut : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  DisplaceInterpolatePut(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, const double *disp1, OutputVoxelType *out)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyDisplacement(disp1);
    this->WorldToInput();
    this->template InterpolatePut(i, j, k, out);
  }

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, const double *o2w, const double *disp1, OutputVoxelType *out)
  {
    this->OutputToWorld(o2w);
    this->ApplyDisplacement(disp1);
    this->WorldToInput();
    this->template InterpolatePut(i, j, k, out);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct DisplaceInterpolatePutAsDouble : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  DisplaceInterpolatePutAsDouble(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  void operator ()(int i, int j, int k, int, const double *disp1, void *)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyDisplacement(disp1);
    this->WorldToInput();
    this->InterpolatePutAsDouble(i, j, k);
  }

  void operator ()(int i, int j, int k, int, const double *o2w, const double *disp1, void *)
  {
    this->OutputToWorld(o2w);
    this->ApplyDisplacement(disp1);
    this->WorldToInput();
    this->InterpolatePutAsDouble(i, j, k);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct DisplaceInterpolateRescalePut : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  DisplaceInterpolateRescalePut(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, const double *disp1, OutputVoxelType *out)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyDisplacement(disp1);
    this->WorldToInput();
    this->template InterpolateRescalePut(i, j, k, out);
  }

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, const double *o2w, const double *disp1, OutputVoxelType *out)
  {
    this->OutputToWorld(o2w);
    this->ApplyDisplacement(disp1);
    this->WorldToInput();
    this->template InterpolateRescalePut(i, j, k, out);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct DisplaceInterpolateRescalePutAsDouble : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  DisplaceInterpolateRescalePutAsDouble(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  void operator ()(int i, int j, int k, int, const double *disp1, void *)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyDisplacement(disp1);
    this->WorldToInput();
    this->InterpolateRescalePutAsDouble(i, j, k);
  }

  void operator ()(int i, int j, int k, int, const double *o2w, const double *disp1, void *)
  {
    this->OutputToWorld(o2w);
    this->ApplyDisplacement(disp1);
    this->WorldToInput();
    this->InterpolateRescalePutAsDouble(i, j, k);
  }
};

// -----------------------------------------------------------------------------
// Using composition of transformations
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct FluidTransformInterpolatePut : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  FluidTransformInterpolatePut(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, OutputVoxelType *out)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyTransformations();
    this->WorldToInput();
    this->template InterpolatePut(i, j, k, out);
  }

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, const double *o2w, OutputVoxelType *out)
  {
    this->OutputToWorld(o2w);
    this->ApplyTransformations();
    this->WorldToInput();
    this->template InterpolatePut(i, j, k, out);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct FluidTransformInterpolatePutAsDouble : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  FluidTransformInterpolatePutAsDouble(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  void operator ()(int i, int j, int k, int, void *)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyTransformations();
    this->WorldToInput();
    this->InterpolatePutAsDouble(i, j, k);
  }

  void operator ()(int i, int j, int k, int, const double *o2w, void *)
  {
    this->OutputToWorld(o2w);
    this->ApplyTransformations();
    this->WorldToInput();
    this->InterpolatePutAsDouble(i, j, k);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct FluidTransformInterpolateRescalePut : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  FluidTransformInterpolateRescalePut(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, OutputVoxelType *out)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyTransformations();
    this->WorldToInput();
    this->template InterpolateRescalePut(i, j, k, out);
  }

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, const double *o2w, OutputVoxelType *out)
  {
    this->OutputToWorld(o2w);
    this->ApplyTransformations();
    this->WorldToInput();
    this->template InterpolateRescalePut(i, j, k, out);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct FluidTransformInterpolateRescalePutAsDouble : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  FluidTransformInterpolateRescalePutAsDouble(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  void operator ()(int i, int j, int k, int, void *)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyTransformations();
    this->WorldToInput();
    this->InterpolateRescalePutAsDouble(i, j, k);
  }

  void operator ()(int i, int j, int k, int, const double *o2w, void *)
  {
    this->OutputToWorld(o2w);
    this->ApplyTransformations();
    this->WorldToInput();
    this->InterpolateRescalePutAsDouble(i, j, k);
  }
};

// -----------------------------------------------------------------------------
// Using composition of displacement and transformation
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct DisplaceTransformInterpolatePut : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  DisplaceTransformInterpolatePut(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {
    if (this->_Transformation1 && !this->_Transformation2) {
      this->_Transformation2 = this->_Transformation1;
    }
  }

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, const double *disp1, OutputVoxelType *out)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyDisplacementAndTransformation(disp1);
    this->WorldToInput();
    this->template InterpolatePut(i, j, k, out);
  }

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, const double *o2w, const double *disp1, OutputVoxelType *out)
  {
    this->OutputToWorld(o2w);
    this->ApplyDisplacementAndTransformation(disp1);
    this->WorldToInput();
    this->template InterpolatePut(i, j, k, out);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct DisplaceTransformInterpolatePutAsDouble : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  DisplaceTransformInterpolatePutAsDouble(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {
    if (this->_Transformation1 && !this->_Transformation2) {
      this->_Transformation2 = this->_Transformation1;
    }
  }

  void operator ()(int i, int j, int k, int, const double *disp1, void *)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyDisplacementAndTransformation(disp1);
    this->WorldToInput();
    this->InterpolatePutAsDouble(i, j, k);
  }

  void operator ()(int i, int j, int k, int, const double *o2w, const double *disp1, void *)
  {
    this->OutputToWorld(o2w);
    this->ApplyDisplacementAndTransformation(disp1);
    this->WorldToInput();
    this->InterpolatePutAsDouble(i, j, k);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct DisplaceTransformInterpolateRescalePut : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  DisplaceTransformInterpolateRescalePut(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {
    if (this->_Transformation1 && !this->_Transformation2) {
      this->_Transformation2 = this->_Transformation1;
    }
  }

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, const double *disp1, OutputVoxelType *out)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyDisplacementAndTransformation(disp1);
    this->WorldToInput();
    this->template InterpolateRescalePut(i, j, k, out);
  }

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, const double *o2w, const double *disp1, OutputVoxelType *out)
  {
    this->OutputToWorld(o2w);
    this->ApplyDisplacementAndTransformation(disp1);
    this->WorldToInput();
    this->template InterpolateRescalePut(i, j, k, out);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct DisplaceTransformInterpolateRescalePutAsDouble : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  DisplaceTransformInterpolateRescalePutAsDouble(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {
    if (this->_Transformation1 && !this->_Transformation2) {
      this->_Transformation2 = this->_Transformation1;
    }
  }

  void operator ()(int i, int j, int k, int, const double *disp1, void *)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyDisplacementAndTransformation(disp1);
    this->WorldToInput();
    this->InterpolateRescalePutAsDouble(i, j, k);
  }

  void operator ()(int i, int j, int k, int, const double *o2w, const double *disp1, void *)
  {
    this->OutputToWorld(o2w);
    this->ApplyDisplacementAndTransformation(disp1);
    this->WorldToInput();
    this->InterpolateRescalePutAsDouble(i, j, k);
  }
};

// -----------------------------------------------------------------------------
// Using additive composition of displacements
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct AddDisplacementsInterpolatePut : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  AddDisplacementsInterpolatePut(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, const double *disp1, const double *disp2, OutputVoxelType *out)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyDisplacements(disp1, disp2);
    this->WorldToInput();
    this->template InterpolatePut(i, j, k, out);
  }

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, const double *o2w, const double *disp1, const double *disp2, OutputVoxelType *out)
  {
    this->OutputToWorld(o2w);
    this->ApplyDisplacements(disp1, disp2);
    this->WorldToInput();
    this->template InterpolatePut(i, j, k, out);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct AddDisplacementsInterpolatePutAsDouble : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  AddDisplacementsInterpolatePutAsDouble(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  void operator ()(int i, int j, int k, int, const double *disp1, const double *disp2, void *)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyDisplacements(disp1, disp2);
    this->WorldToInput();
    this->InterpolatePutAsDouble(i, j, k);
  }

  void operator ()(int i, int j, int k, int, const double *o2w, const double *disp1, const double *disp2, void *)
  {
    this->OutputToWorld(o2w);
    this->ApplyDisplacements(disp1, disp2);
    this->WorldToInput();
    this->InterpolatePutAsDouble(i, j, k);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct AddDisplacementsInterpolateRescalePut : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  AddDisplacementsInterpolateRescalePut(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, const double *disp1, const double *disp2, OutputVoxelType *out)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyDisplacements(disp1, disp2);
    this->WorldToInput();
    this->template InterpolateRescalePut(i, j, k, out);
  }

  template <class OutputVoxelType>
  void operator ()(int i, int j, int k, int, const double *o2w, const double *disp1, const double *disp2, OutputVoxelType *out)
  {
    this->OutputToWorld(o2w);
    this->ApplyDisplacements(disp1, disp2);
    this->WorldToInput();
    this->template InterpolateRescalePut(i, j, k, out);
  }
};

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction = InterpolateImageFunction, class InputDomain = InterpolationDomain::Foreground>
struct AddDisplacementsInterpolateRescalePutAsDouble : public BaseTransform<InterpolateImageFunction, InputDomain>
{
  AddDisplacementsInterpolateRescalePutAsDouble(const TransformData &data)
    : BaseTransform<InterpolateImageFunction, InputDomain>(data)
  {}

  void operator ()(int i, int j, int k, int, const double *disp1, const double *disp2, void *)
  {
    this->OutputToWorld(i, j, k);
    this->ApplyDisplacements(disp1, disp2);
    this->WorldToInput();
    this->InterpolateRescalePutAsDouble(i, j, k);
  }

  void operator ()(int i, int j, int k, int, const double *o2w, const double *disp1, const double *disp2, void *)
  {
    this->OutputToWorld(o2w);
    this->ApplyDisplacements(disp1, disp2);
    this->WorldToInput();
    this->InterpolateRescalePutAsDouble(i, j, k);
  }
};


} } // namespace mirtk::MultipleVoxelTransformation

#endif // MIRTK_MultipleVoxelTransformation_H
