/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2010-2015 Imperial College London
 * Copyright 2010      Wenzhe Shi
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
 
#include "mirtk/ShapeBasedInterpolateImageFunction.h"

#include "mirtk/Math.h"
#include "mirtk/Resampling.h"
#include "mirtk/EuclideanDistanceTransform.h"
#include "mirtk/LinearInterpolateImageFunction.hxx"

#include "mirtk/CommonExport.h"


namespace mirtk {


// Global "verbose" flag (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int verbose;


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
ShapeBasedInterpolateImageFunction::ShapeBasedInterpolateImageFunction()
{
}

// -----------------------------------------------------------------------------
ShapeBasedInterpolateImageFunction::~ShapeBasedInterpolateImageFunction()
{
}

// -----------------------------------------------------------------------------
void ShapeBasedInterpolateImageFunction::Refine()
{
  typedef EuclideanDistanceTransform<RealPixel> DistanceTransformType;

  int    i, j, k, l;
  double min, max, current, icurrent, dcurrent, rcurrent, sum, sumcount;

  // Initialization
  sum = .0;

  // For every intensity value
  Input()->GetMinMaxAsDouble(&min, &max);

  // Initialize _rcdmap
  for (l = 0; l < _rcdmap.T(); ++l)
  for (k = 0; k < _rcdmap.Z(); ++k)
  for (j = 0; j < _rcdmap.Y(); ++j)
  for (i = 0; i < _rcdmap.X(); ++i) {
    _rcdmap.PutAsDouble(i, j, k, l, 30);
  }

  for (current = min; current <= max; current++) {
    // Threshold
    sumcount = .0;
    for (l = 0; l < _tinput.T(); ++l)
    for (k = 0; k < _tinput.Z(); ++k)
    for (j = 0; j < _tinput.Y(); ++j)
    for (i = 0; i < _tinput.X(); ++i) {
      icurrent = Input()->GetAsDouble(i, j, k, l);
      if (icurrent < current || icurrent >= (current + 1)) {
        _tinput(i, j, k, l) = 0;
      } else {
        _tinput(i, j, k, l) = 1;
        ++sumcount;
      }
    }

    // Calculate EDT
    if (verbose && iround(current) % 20 == 0) {
      if (max < current + 19) {
        cout << "Doing outside DT for value == : "<< current << " to " << max << endl;
        cout << "Doing inside  DT for value == : "<< current << " to " << max << endl;
      } else {
        cout << "Doing outside DT for value == : "<< current << " to " << current + 19.0 << endl;
        cout << "Doing inside  DT for value == : "<< current << " to " << current + 19.0 << endl;
      }
    }

    if (sumcount > .0) {
      // Dmap _tinput to _dmap
      {
        RealImage inputA, inputB, outputA, outputB;

        // Default mode
        DistanceTransformType edt(DistanceTransformType::DT_3D);

        // Threshold image
        inputA = _tinput;
        inputB = _tinput;
        for (l = 0; l < _tinput.T(); ++l)
        for (k = 0; k < _tinput.Z(); ++k)
        for (j = 0; j < _tinput.Y(); ++j)
        for (i = 0; i < _tinput.X(); ++i) {
          if (_tinput(i, j, k, l) > .5) {
            inputA(i, j, k, l) = 1;
            inputB(i, j, k, l) = 0;
          } else {
            inputA(i, j, k, l) = 0;
            inputB(i, j, k, l) = 1;
          }
        }

        edt.Input (& inputA);
        edt.Output(&outputA);
        edt.Run();		  
        edt.Input (& inputB);
        edt.Output(&outputB);
        edt.Run();

        for (l = 0; l < _tinput.T(); ++l)
        for (k = 0; k < _tinput.Z(); ++k)
        for (j = 0; j < _tinput.Y(); ++j)
        for (i = 0; i < _tinput.X(); ++i) {
          _dmap(i, j, k, l) = sqrt(outputA(i, j, k, l)) - sqrt(outputB(i, j, k, l));
        }
      }

      // Linear Interpolate Dmap _dmap to _rdmap
      {
        double x, y, z, t;
        LinearInterpolateImageFunction interpolator;
        interpolator.Input(&_dmap);
        interpolator.Initialize();

        for (l = 0; l < _rdmap.T(); ++l)
        for (k = 0; k < _rdmap.Z(); ++k)
        for (j = 0; j < _rdmap.Y(); ++j)
        for (i = 0; i < _rdmap.X(); ++i) {
          x = i, y = j, z = k, t = l;
          _rdmap.ImageToWorld(x, y, z);
          _dmap .WorldToImage(x, y, z);
          _rdmap.PutAsDouble(i, j, k, l, interpolator.Evaluate(x, y, z, t));
        }
      }

      // Put value back to Resampled Image _rinput < 0 if _rinput == 0 let the neighbor vote.
      {
        for (l = 0; l < _rdmap.T(); ++l)
        for (k = 0; k < _rdmap.Z(); ++k)
        for (j = 0; j < _rdmap.Y(); ++j)
        for (i = 0; i < _rdmap.X(); ++i) {
          icurrent = _rinput.GetAsDouble(i, j, k, l);
          dcurrent = _rdmap .GetAsDouble(i, j, k, l);
          rcurrent = _rcdmap.GetAsDouble(i, j, k, l);
          if (dcurrent < 0 && current >= icurrent && dcurrent < rcurrent) {
            _rcdmap.PutAsDouble(i, j, k, l, dcurrent);
          } else if (dcurrent <= .0 && current >= icurrent && dcurrent < rcurrent) {
            sum = .0, sumcount = .0;
            if (i > 0) {
              sum += _rdmap.GetAsDouble(i-1, j, k, l);
              ++sumcount;
            } else if (i < _rinput.X() - 1) {
              sum += _rdmap.GetAsDouble(i+1, j, k, l);
              ++sumcount;
            }
            if (j > 0) {
              sum += _rdmap.GetAsDouble(i, j-1, k, l);
              ++sumcount;
            } else if (j < _rinput.Y() - 1) {
              sum += _rdmap.GetAsDouble(i, j+1, k, l);
              ++sumcount;
            }
            if (k > 0) {
              sum += _rdmap.GetAsDouble(i, j, k-1, l);
              ++sumcount;
            } else if (k < _rinput.Z() - 1) {
              sum += _rdmap.GetAsDouble(i, j, k+1, l);
              ++sumcount;
            }
            sum = sum/sumcount;
            if (sum <= .0) {
              _rcdmap.PutAsDouble(i, j, k, l, dcurrent);
            }
          } else if (dcurrent > .0 && icurrent <= current && dcurrent < rcurrent) {
            _rinput.PutAsDouble(i, j, k, l, current);
            _rcdmap.PutAsDouble(i, j, k, l, dcurrent);
          } else if (dcurrent >= .0 && icurrent <= current && dcurrent < rcurrent) {
            if (sum > .0) {
              _rinput.PutAsDouble(i, j, k, l, current);
              _rcdmap.PutAsDouble(i, j, k, l, dcurrent);
            }
          }
        }
      }
    }
  } // End for
}

// ----------------------------------------------------------------------------
void ShapeBasedInterpolateImageFunction::Initialize(bool coeff)
{
  typedef EuclideanDistanceTransform<RealPixel> DistanceTransformType;

  // Initialize base class
  InterpolateImageFunction::Initialize(coeff);

  double xsize, ysize, zsize, size;
  int    new_x, new_y, new_z, i, j, k, l, labelcount;
  double xaxis[3], yaxis[3], zaxis[3];
  double new_xsize, new_ysize, new_zsize;
  double old_xsize, old_ysize, old_zsize;
  double min, max, current, sum, sumcount;

  // Determine minimum voxel size
  Input()->GetPixelSize(&xsize, &ysize, &zsize);
  size = xsize;
  size = (size < ysize) ? size : ysize;
  size = (size < zsize) ? size : zsize;
  if (size > 1) size = 1;
  cerr << "Create images with isotropic voxel size (in mm): "<< size << endl;

  // Create _rinput _rdmap
  _dmap   = RealImage(Input()->Attributes());
  _tinput = RealImage(Input()->Attributes());

  // Determine the old dimensions of the image
  Input()->GetPixelSize(&old_xsize, &old_ysize, &old_zsize);

  // Determine the new dimensions of the image
  new_x = int(Input()->X() * old_xsize / size);
  new_y = int(Input()->Y() * old_ysize / size);
  new_z = int(Input()->Z() * old_zsize / size);

  // Determine the new voxel dimensions
  if (new_x < 1) {
    new_x     =  1;
    new_xsize =  old_xsize;
  } else {
    new_xsize = size;
  }
  if (new_y < 1) {
    new_y     =  1;
    new_ysize =  old_ysize;
  } else {
    new_ysize = size;
  }
  if (new_z < 1) {
    new_z     =  1;
    new_zsize =  old_zsize;
  } else {
    new_zsize = size;
  }

  // Allocate new image
  _rinput = RealImage(new_x, new_y, new_z, Input()->T());
  _rdmap  = RealImage(new_x, new_y, new_z, Input()->T());
  _rcdmap = RealImage(new_x, new_y, new_z, Input()->T());

  // Set new voxel size
  _rinput.PutPixelSize(new_xsize, new_ysize, new_zsize);
  _rdmap .PutPixelSize(new_xsize, new_ysize, new_zsize);
  _rcdmap.PutPixelSize(new_xsize, new_ysize, new_zsize);

  // Set new orientation
  Input()->GetOrientation(xaxis, yaxis, zaxis);
  _rinput.PutOrientation(xaxis, yaxis, zaxis);
  _rdmap .PutOrientation(xaxis, yaxis, zaxis);
  _rcdmap.PutOrientation(xaxis, yaxis, zaxis);

  // Set new origin
  _rinput.PutOrigin(Input()->GetOrigin());
  _rdmap .PutOrigin(Input()->GetOrigin());
  _rcdmap.PutOrigin(Input()->GetOrigin());

  // For every intensity value
  Input()->GetMinMaxAsDouble(&min, &max);

  labelcount = 0;
  for (current = min; current <= max; ++current) {
    // Threshold
    sumcount = .0;
    for (l = 0; l < _tinput.T(); ++l)
    for (k = 0; k < _tinput.Z(); ++k)
    for (j = 0; j < _tinput.Y(); ++j)
    for (i = 0; i < _tinput.X(); ++i) {
      if (Input()->GetAsDouble(i, j, k, l) < current /* || Input()->GetAsDouble(i, j, k, l) > current */) {
        _tinput(i, j, k, l) = 0;
      } else {
        _tinput(i, j, k, l) = 1;
      }
      if (Input()->GetAsDouble(i, j, k, l) >= current &&
          Input()->GetAsDouble(i, j, k, l) <  current + 1.0) {
        ++sumcount;
      }
    }

    if (verbose && iround(current) % 20 == 0) {
      if (max < current + 19.0) {
        cout << "Doing outside DT for value >= : "<< current << " to " << max << endl;
        cout << "Doing inside  DT for value >= : "<< current << " to " << max << endl;
      } else {
        cout << "Doing outside DT for value >= : "<< current << " to " << current + 19.0 << endl;
        cout << "Doing inside  DT for value >= : "<< current << " to " << current + 19.0 << endl;
      }
    }

    if (sumcount > .0) {
      ++labelcount;
      // Dmap _tinput to _dmap
      {
        RealImage inputA, inputB, outputA, outputB;

        // Default mode
        DistanceTransformType edt(DistanceTransformType::DT_3D);

        // Threshold image
        inputA = _tinput;
        inputB = _tinput;
        for (l = 0; l < _tinput.T(); ++l)
        for (k = 0; k < _tinput.Z(); ++k)
        for (j = 0; j < _tinput.Y(); ++j)
        for (i = 0; i < _tinput.X(); ++i) {
          if (_tinput(i, j, k, l) > .5) {
            inputA(i, j, k, l) = 1;
            inputB(i, j, k, l) = 0;
          } else {
            inputA(i, j, k, l) = 0;
            inputB(i, j, k, l) = 1;
          }
        }

        // Calculate EDT
        edt.Input (& inputA);
        edt.Output(&outputA);
        edt.Run();		  
        edt.Input (& inputB);
        edt.Output(&outputB);
        edt.Run();

        for (l = 0; l < _tinput.T(); ++l)
        for (k = 0; k < _tinput.Z(); ++k)
        for (j = 0; j < _tinput.Y(); ++j)
        for (i = 0; i < _tinput.X(); ++i) {
          _dmap(i, j, k, l)  = sqrt(outputA(i, j, k, l)) - sqrt(outputB(i, j, k, l));
        }
      }

      // Linear Interpolate Dmap _dmap to _rdmap
      {
        double x, y, z, t;
        LinearInterpolateImageFunction interpolator;
        interpolator.Input(&_dmap);
        interpolator.Initialize();

        for (l = 0; l < _rdmap.T(); ++l)
        for (k = 0; k < _rdmap.Z(); ++k)
        for (j = 0; j < _rdmap.Y(); ++j)
        for (i = 0; i < _rdmap.X(); ++i) {
          x = i, y = j, z = k, t = l;
          _rdmap.ImageToWorld(x, y, z);
          _dmap .WorldToImage(x, y, z);
          _rdmap.PutAsDouble(i, j, k, l, interpolator.Evaluate(x, y, z, t));
        }
      }

      // Put value back to Resampled Image _rinput < 0 if _rinput == 0 let the neighbor vote.
      {
        for (l = 0; l < _rdmap.T(); ++l)
        for (k = 0; k < _rdmap.Z(); ++k)
        for (j = 0; j < _rdmap.Y(); ++j)
        for (i = 0; i < _rdmap.X(); ++i) {
          if (_rdmap.GetAsDouble(i, j, k, l) < .0) {
            _rinput.PutAsDouble(i, j, k, l, current);
          } else if (_rdmap.GetAsDouble(i, j, k, l) <= .0) {
            sum = .0, sumcount = .0;
            if (i > 0) {
              sum += _rdmap.GetAsDouble(i-1, j, k, l);
              ++sumcount;
            } else if (i < _rinput.X() - 1) {
              sum += _rdmap.GetAsDouble(i+1, j, k, l);
              ++sumcount;
            }
            if (j > 0) {
              sum += _rdmap.GetAsDouble(i, j-1, k, l);
              ++sumcount;
            } else if (j < _rinput.Y() - 1) {
              sum += _rdmap.GetAsDouble(i, j+1, k, l);
              ++sumcount;
            }
            if (k > 0) {
              sum += _rdmap.GetAsDouble(i, j, k-1, l);
              ++sumcount;
            } else if (k < _rinput.Z() - 1) {
              sum += _rdmap.GetAsDouble(i, j, k+1, l);
              ++sumcount;
            }
            sum = sum/sumcount;
            if (sum <= .0) {
              _rinput.PutAsDouble(i, j, k, l, current);
            }
          }
        }
      }
    }

  } // End for

  // Refine to fix the union property
  if (labelcount > 3 && labelcount < 50) Refine();

  // Instantiate internal interpolator
  BaseImage *input = &_rinput;

  _nn_interpolator.Input(input);
  _nn_interpolator.Initialize();

  _linear_interpolator.Input(input);
  _linear_interpolator.Initialize();

  // Domain on which interpolation is defined
  _nn_interpolator.Inside(this->_x1, this->_y1, this->_z1,
                          this->_x2, this->_y2, this->_z2);
  _rinput .ImageToWorld(this->_x1, this->_y1, this->_z1);
  Input()->WorldToImage(this->_x1, this->_y1, this->_z1);
  _rinput .ImageToWorld(this->_x2, this->_y2, this->_z2);
  Input()->WorldToImage(this->_x2, this->_y2, this->_z2);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double ShapeBasedInterpolateImageFunction::EvaluateInside(double x, double y, double z, double t) const
{
  Input()->ImageToWorld(x, y, z);
  _rinput.WorldToImage(x, y, z);
  return _nn_interpolator.EvaluateInside(x, y, z, t);
}

// -----------------------------------------------------------------------------
double ShapeBasedInterpolateImageFunction::EvaluateOutside(double x, double y, double z, double t) const
{
  Input()->ImageToWorld(x, y, z);
  _rinput.WorldToImage(x, y, z);
  return _nn_interpolator.EvaluateOutside(x, y, z, t);
}

// -----------------------------------------------------------------------------
double ShapeBasedInterpolateImageFunction::EvaluateWithPaddingInside(double x, double y, double z, double t) const
{
  Input()->ImageToWorld(x, y, z);
  _rinput.WorldToImage(x, y, z);
  return _nn_interpolator.EvaluateWithPaddingInside(x, y, z, t);
}

// -----------------------------------------------------------------------------
double ShapeBasedInterpolateImageFunction::EvaluateWithPaddingOutside(double x, double y, double z, double t) const
{
  Input()->ImageToWorld(x, y, z);
  _rinput.WorldToImage(x, y, z);
  return _nn_interpolator.EvaluateWithPaddingOutside(x, y, z, t);
}

// -----------------------------------------------------------------------------
double ShapeBasedInterpolateImageFunction::EvaluateLinear(double x, double y, double z, double t) const
{
  if (IsInside(x, y, z)) return this->EvaluateInsideLinear (x, y, z, t);
  else                   return this->EvaluateOutsideLinear(x, y, z, t);
}

// -----------------------------------------------------------------------------
double ShapeBasedInterpolateImageFunction::EvaluateInsideLinear(double x, double y, double z, double t) const
{
  Input()->ImageToWorld(x, y, z);
  _rinput.WorldToImage(x, y, z);
  return _linear_interpolator.EvaluateInside(x, y, z, t);
}

// -----------------------------------------------------------------------------
double ShapeBasedInterpolateImageFunction::EvaluateOutsideLinear(double x, double y, double z, double t) const
{
  Input()->ImageToWorld(x, y, z);
  _rinput.WorldToImage(x, y, z);
  return _linear_interpolator.EvaluateOutside(x, y, z, t);
}


} // namespace mirtk
