/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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

#include "mirtk/ImageTransformation.h"

#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"

#include "mirtk/HomogeneousTransformation.h"
#include "mirtk/HomogeneousTransformationIterator.h"


namespace mirtk {


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace ImageTransformationUtils {


// -----------------------------------------------------------------------------
/// Applies any transformation to resample a given source image on a target grid
struct ApplyTransformation
{
  // ---------------------------------------------------------------------------
  // Attributes
  const BaseImage                *_Input;
  const ImageFunction            *_Interpolator;
  const Transformation           *_Transformation;
  const InterpolateImageFunction *_DisplacementField;
  BaseImage                      *_Output;

  bool   _Invert;
  bool   _TwoD;
  double _InputMin;
  double _InputMax;
  double _ScaleFactor;
  double _Offset;
  double _TargetPaddingValue;
  double _SourcePaddingValue;
  double _InputTimeOffset;
  double _OutputTimeOffset;

  mutable int    _OutputFrame;
  mutable double _OutputTime;
  mutable int    _InputFrame;
  mutable double _InputTime;

  int _NumberOfSingularPoints;

  // ---------------------------------------------------------------------------
  /// Default constructor
  ApplyTransformation()
  :
    _Input             (NULL),
    _Interpolator      (NULL),
    _Transformation    (NULL),
    _DisplacementField (NULL),
    _Output            (NULL),
    _Invert            (false),
    _TwoD              (false),
    _InputMin          (.0),
    _InputMax          (-1),
    _ScaleFactor       (1),
    _Offset            (0),
    _TargetPaddingValue(-1),
    _SourcePaddingValue(-1),
    _InputTimeOffset   (0),
    _OutputTimeOffset  (0),
    _OutputFrame       (0),
    _OutputTime        (-1),
    _InputFrame        (0),
    _InputTime         (0),
    _NumberOfSingularPoints(0)
  {}

  // ---------------------------------------------------------------------------
  /// Split constructor
  ApplyTransformation(const ApplyTransformation &other, split)
  :
    _Input             (other._Input),
    _Interpolator      (other._Interpolator),
    _Transformation    (other._Transformation),
    _DisplacementField (other._DisplacementField),
    _Output            (other._Output),
    _Invert            (other._Invert),
    _TwoD              (other._TwoD),
    _InputMin          (other._InputMin),
    _InputMax          (other._InputMax),
    _ScaleFactor       (other._ScaleFactor),
    _Offset            (other._Offset),
    _TargetPaddingValue(other._TargetPaddingValue),
    _SourcePaddingValue(other._SourcePaddingValue),
    _InputTimeOffset   (other._InputTimeOffset),
    _OutputTimeOffset  (other._OutputTimeOffset),
    _OutputFrame       (other._OutputFrame),
    _OutputTime        (other._OutputTime),
    _InputFrame        (other._InputFrame),
    _InputTime         (other._InputTime),
    _NumberOfSingularPoints(0)
  {}

  // ---------------------------------------------------------------------------
  void join(const ApplyTransformation &other)
  {
    _NumberOfSingularPoints += other._NumberOfSingularPoints;
  }

  // ---------------------------------------------------------------------------
  double Clamp(double v) const
  {
    if (_InputMin <= _InputMax) {
      if      (v < _InputMin) return _InputMin;
      else if (v > _InputMax) return _InputMax;
    }
    return v;
  }

  // ---------------------------------------------------------------------------
  /// Applies the given transformation within given output region for frame _T
  void operator() (const blocked_range3d<int> &r)
  {
    double value, x, y, z, u, v, w, disp[3] = {.0, .0, .0};

    for (int k = r.pages().begin(); k != r.pages().end(); ++k)
    for (int j = r.rows ().begin(); j != r.rows ().end(); ++j)
    for (int i = r.cols ().begin(); i != r.cols ().end(); ++i) {
      if (_Output->GetAsDouble(i, j, k, _OutputFrame) > _TargetPaddingValue) {
        // Transform point into world coordinates
        x = i, y = j, z = k;
        _Output->ImageToWorld(x, y, z);
        // Transform point
        if (_DisplacementField) {
          u = x, v = y, w = z;
          _DisplacementField->WorldToImage  (u, v, w);
          _DisplacementField->Evaluate(disp, u, v, w);
          x += disp[0], y += disp[1], z += disp[2];
        } else {
          if (_Invert) {
            if (!_Transformation->Inverse(x, y, z, _InputTime, _OutputTime)) {
              ++_NumberOfSingularPoints;
            }
          } else {
            _Transformation->Transform(x, y, z, _InputTime, _OutputTime);
          }
        }
        // Transform point into image coordinates
        _Input->WorldToImage(x, y, z);
        // Check whether transformed point is in FOV of input
        if (-0.5 < x && x < static_cast<double>(_Input->X()) - 0.5 &&
            -0.5 < y && y < static_cast<double>(_Input->Y()) - 0.5) {
          if (_TwoD) {
            value = _Interpolator->Evaluate(x, y, k, _InputFrame);
            value = _ScaleFactor * Clamp(value) + _Offset;
          } else if (-0.5 < z && z < static_cast<double>(_Input->Z()) - 0.5) {
            value = _Interpolator->Evaluate(x, y, z, _InputFrame);
            value = _ScaleFactor * Clamp(value) + _Offset;
          } else {
            value = _SourcePaddingValue;
          }
        } else {
          value = _SourcePaddingValue;
        }
      } else {
        value = _SourcePaddingValue;
      }
      _Output->PutAsDouble(i, j, k, _OutputFrame, value);
    }
  }

  // ---------------------------------------------------------------------------
  /// Applies the given transformation for the given range of output frames
  void operator() ()
  {
    _NumberOfSingularPoints = 0;
    for (_OutputFrame = 0; _OutputFrame < _Output->T(); ++_OutputFrame) {
      _OutputTime = _Output->ImageToTime(_OutputFrame);
      _InputFrame = ((_Input->T() > 1) ? iround(_Input->TimeToImage(_OutputTime)) : 0);
      _InputTime  = _Input->ImageToTime(_InputFrame);

      _OutputTime += _OutputTimeOffset;
      _InputTime  += _InputTimeOffset;

      if (0 <= _InputFrame && _InputFrame < _Input->T()) {
        blocked_range3d<int> voxels(0, _Output->Z(),
                                    0, _Output->Y(),
                                    0, _Output->X());
        ApplyTransformation body(*this);
        parallel_reduce(voxels, body);
        _NumberOfSingularPoints += body._NumberOfSingularPoints;
      } else {
        for (int k = 0; k < _Output->Z(); ++k)
        for (int j = 0; j < _Output->Y(); ++j)
        for (int i = 0; i < _Output->X(); ++i) {
          _Output->PutAsDouble(i, j, k, _OutputFrame, _SourcePaddingValue);
        }
      }
    }
  }

};

// -----------------------------------------------------------------------------
/// Applies linear homogeneous coordinate transformation
struct ApplyHomogeneousTransformation
{
  const BaseImage                 *_Input;
  const ImageFunction             *_Interpolator;
  const HomogeneousTransformation *_Transformation;
  BaseImage                       *_Output;

  bool   _Invert;
  bool   _TwoD;
  double _InputMin;
  double _InputMax;
  double _ScaleFactor;
  double _Offset;
  double _TargetPaddingValue;
  double _SourcePaddingValue;
  int    _InputFrame;
  int    _OutputFrame;

  // ---------------------------------------------------------------------------
  double Clamp(double v) const
  {
    if (_InputMin <= _InputMax) {
      if      (v < _InputMin) return _InputMin;
      else if (v > _InputMax) return _InputMax;
    }
    return v;
  }

  // ---------------------------------------------------------------------------
  /// Applies the given transformation within given output region for frame _T
  void operator()(const blocked_range<int> &r) const
  {
    double value;
    HomogeneousTransformationIterator it(_Transformation);
    it.Initialize(_Output, _Input, .0, .0, static_cast<double>(r.begin()), _Invert);
    for (int k = r.begin(); k != r.end();     ++k, it.NextZ())
    for (int j = 0;         j < _Output->Y(); ++j, it.NextY())
    for (int i = 0;         i < _Output->X(); ++i, it.NextX()) {
      if (_Output->GetAsDouble(i, j, k, _OutputFrame) > _TargetPaddingValue) {
        if (-.5 < it._x && it._x < _Input->X() - .5 &&
            -.5 < it._y && it._y < _Input->Y() - .5) {
          if (_TwoD) {
            value = _Interpolator->Evaluate(it._x, it._y, k, _InputFrame);
            value = _ScaleFactor * Clamp(value) + _Offset;
          } else if (-.5 < it._z && it._z < _Input->Z() - .5) {
            value = _Interpolator->Evaluate(it._x, it._y, it._z, _InputFrame);
            value = _ScaleFactor * Clamp(value) + _Offset;
          } else {
            value = _SourcePaddingValue;
          }
        } else {
          value = _SourcePaddingValue;
        }
      } else {
      	value = _SourcePaddingValue;
      }
      _Output->PutAsDouble(i, j, k, _OutputFrame, value);
    }
  }

  // ---------------------------------------------------------------------------
  /// Applies the given transformation for the given range of output frames
  void operator() ()
  {
    for (_OutputFrame = 0; _OutputFrame < _Output->T(); ++_OutputFrame) {
      double time = _Output->ImageToTime(_OutputFrame);
      _InputFrame = ((_Input->T() > 1) ? iround(_Input->TimeToImage(time)) : 0);
      if (0 <= _InputFrame && _InputFrame < _Input->T()) {
        ApplyHomogeneousTransformation body(*this);
        parallel_for(blocked_range<int>(0, _Output->Z()), body);
      } else {
        for (int k = 0; k < _Output->Z(); ++k)
        for (int j = 0; j < _Output->Y(); ++j)
        for (int i = 0; i < _Output->X(); ++i) {
          _Output->PutAsDouble(i, j, k, _OutputFrame, _SourcePaddingValue);
        }
      }
    }
  }
};


} // namespace ImageTransformationUtils
using namespace ImageTransformationUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
ImageTransformation::ImageTransformation()
:
  _Input(NULL),
  _Output(NULL),
  _Transformation(NULL),
  _Interpolator(NULL),
  _TargetPaddingValue(voxel_limits<double>::min()),
  _SourcePaddingValue(.0),
  _ScaleFactor(1.0),
  _Offset(.0),
  _InputTimeOffset(.0),
  _OutputTimeOffset(.0),
  _Invert(false),
  _TwoD(false),
  _Cache(NULL),
  _CacheOwner(false),
  _CacheInterpolation(Interpolation_Linear),
  _CacheExtrapolation(Extrapolation_Const),
  _DisplacementField(NULL),
  _NumberOfSingularPoints(0)
{
}

// -----------------------------------------------------------------------------
ImageTransformation::~ImageTransformation()
{
  if (_CacheOwner) Delete(_Cache);
  Delete(_DisplacementField);
}

// =============================================================================
// Setters
// =============================================================================

// -----------------------------------------------------------------------------
void ImageTransformation::Transformation(const class Transformation *transformation)
{
  if (_Transformation != transformation) {
    _Transformation = transformation;
    if (_Cache) _Cache->Modified(true);
  }
}

// -----------------------------------------------------------------------------
void ImageTransformation::Cache(ImageTransformationCache *cache)
{
  if (_Cache != cache) {
    if (_CacheOwner) Delete(_Cache);
    Delete(_DisplacementField);
    _Cache      = cache;
    _CacheOwner = false;
  }
}

// -----------------------------------------------------------------------------
void ImageTransformation::Invert(bool inv)
{
  if (_Invert != inv)
  {
    if (_Cache) _Cache->Modified(true);
    _Invert = inv;
  }
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void ImageTransformation::Initialize()
{
  // Clean up previous run
  Delete(_DisplacementField);

  // Check inputs and outputs
  if (!_Input) {
    cerr << "ImageTransformation::Initialize: Filter has no input" << endl;
    exit(1);
  }
  if (!_Output) {
    cerr << "ImageTransformation::Initialize: Filter has no output" << endl;
    exit(1);
  }
  if (!_Transformation && (!_Cache || _Cache->IsEmpty() || _Cache->Modified())) {
    cerr << "ImageTransformation::Initialize: Filter has no transformation" << endl;
    exit(1);
  }
  if (!_Interpolator) {
    cerr << "ImageTransformation::Initialize: Filter has no interpolator" << endl;
    exit(1);
  }
  if (_Input->IsEmpty()) {
    cerr << "ImageTransformation::Initialize: Input is empty" << endl;
    exit(1);
  }
  if (_Input == _Output) {
    cerr << "ImageTransformation::Initialize: Input equals output" << endl;
    exit(1);
  }

  _NumberOfSingularPoints = 0;

  if (!_Cache && _Transformation->RequiresCachingOfDisplacements()) {
    _CacheOwner = true;
    _Cache      = new ImageTransformationCache();
    _Cache->Initialize(_Output->Attributes(), 3);
  }

  if (_Cache && !_Cache->IsEmpty()) {

    // Compute and cache displacements
    if (_Cache->Modified()) {
      // TODO: Cache displacements for multiple time intervals.
      const double t0 = _Cache->GetTOrigin();
      const double t  = _Input->GetTOrigin();

      _Cache->Initialize();
      if (_Invert) {
        _NumberOfSingularPoints = _Transformation->InverseDisplacement(*_Cache, t, t0);
      } else {
        _Transformation->Displacement(*_Cache, t, t0);
      }

      _Cache->Modified(false);
    }

    // Initialize cache interpolator
    _DisplacementField = InterpolateImageFunction::New(_CacheInterpolation,
                                                       _CacheExtrapolation,
                                                       _Cache);
    _DisplacementField->Input(_Cache);
    _DisplacementField->Initialize();
  }

  // Initialize input interpolator
  _Interpolator->Input(_Input);
  _Interpolator->Initialize();
}

// -----------------------------------------------------------------------------
void ImageTransformation::Run()
{
  MIRTK_START_TIMING();

  this->Initialize();

  // Range of output intensities should be the same as the one of the input
  double min_intensity, max_intensity;
  _Input->GetMinMaxAsDouble(min_intensity, max_intensity);

  const HomogeneousTransformation *lin;
  if ((lin = dynamic_cast<const HomogeneousTransformation *>(_Transformation))) {
    ApplyHomogeneousTransformation run;
    run._Input              = _Input;
    run._Interpolator       = _Interpolator;
    run._Transformation     = lin;
    run._Output             = _Output;
    run._Invert             = _Invert;
    run._TwoD               = _TwoD;
    run._ScaleFactor        = _ScaleFactor;
    run._Offset             = _Offset;
    run._InputMin           = min_intensity;
    run._InputMax           = max_intensity;
    run._TargetPaddingValue = _TargetPaddingValue;
    run._SourcePaddingValue = _SourcePaddingValue;
    run();
  } else {
    ApplyTransformation run;
    run._Input              = _Input;
    run._Interpolator       = _Interpolator;
    run._Transformation     = _Transformation;
    run._DisplacementField  = _DisplacementField;
    run._Output             = _Output;
    run._Invert             = _Invert;
    run._TwoD               = _TwoD;
    run._ScaleFactor        = _ScaleFactor;
    run._Offset             = _Offset;
    run._InputMin           = min_intensity;
    run._InputMax           = max_intensity;
    run._TargetPaddingValue = _TargetPaddingValue;
    run._SourcePaddingValue = _SourcePaddingValue;
    run._InputTimeOffset    = _InputTimeOffset;
    run._OutputTimeOffset   = _OutputTimeOffset;
    run();
    if (_Invert && !_DisplacementField) {
      _NumberOfSingularPoints = run._NumberOfSingularPoints;
    }
  }

  MIRTK_DEBUG_TIMING(1, "ImageTransformation::Run");
}


} // namespace mirtk
