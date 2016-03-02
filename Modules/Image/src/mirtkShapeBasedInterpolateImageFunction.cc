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
 
#include <mirtkShapeBasedInterpolateImageFunction.h>

#include <mirtkMath.h>
#include <mirtkResampling.h>
#include <mirtkEuclideanDistanceTransform.h>
#include <mirtkLinearInterpolateImageFunction.hxx>


namespace mirtk {


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

  int    x, y, z, t;
  double min, max, current, sum, sumcount, dcurrent;

  // Initialization
  sum = 0;

  // For every intensity value
  Input()->GetMinMaxAsDouble(&min, &max);

  // Initialize _rcdmap
  for (t = 0 ; t < _rcdmap.GetT(); t++) {
    for (z = 0; z < _rcdmap.GetZ(); z++) {
      for (y = 0; y < _rcdmap.GetY(); y++) {
        for (x = 0; x < _rcdmap.GetX(); x++) {
          _rcdmap.PutAsDouble(x,y,z,t,30);
        }
      }
    }
  }

  for (current = min; current <= max; current++) {
    // Threshold
    sumcount = 0;
    for (t = 0; t < _tinput.GetT(); t++) {
      for (z = 0; z < _tinput.GetZ(); z++) {
        for (y = 0; y < _tinput.GetY(); y++) {
          for (x = 0; x < _tinput.GetX(); x++) {
            if ((Input()->GetAsDouble(x, y, z, t) < current) || (Input()->GetAsDouble(x, y, z, t) >= (current + 1))) {
              _tinput(x, y, z, t) = 0;
            } else {
              _tinput(x, y, z, t) = 1;
              sumcount ++;
            }
          }
        }
      }
    }

    // Calculate EDT
    if (iround(current) % 20 == 0) {
      if (max < current + 19) {
        cout << "Doing outside DT for value == : "<< current << " to " << max << endl;
        cout << "Doing inside DT for value == : "<< current << " to " << max << endl;
      } else {
        cout << "Doing outside DT for value == : "<< current << " to " << current+19 << endl;
        cout << "Doing inside DT for value == : "<< current << " to " << current+19 << endl;
      }
    }

    if (sumcount > 0) {
      // Dmap _tinput to _dmap
      {
        RealImage inputA, inputB, outputA, outputB;

        // Default mode
        DistanceTransformType edt(DistanceTransformType::DT_3D);

        // Threshold image
        inputA = _tinput;
        inputB = _tinput;
        for (t = 0; t < _tinput.GetT(); t++) {
          for (z = 0; z < _tinput.GetZ(); z++) {
            for (y = 0; y < _tinput.GetY(); y++) {
              for (x = 0; x < _tinput.GetX(); x++) {
                if (_tinput(x, y, z, t) > 0.5) {
                  inputA(x, y, z, t) = 1;
                  inputB(x, y, z, t) = 0;
                } else {
                  inputA(x, y, z, t) = 0;
                  inputB(x, y, z, t) = 1;
                }
              }
            }
          }
        }

        edt.Input (& inputA);
        edt.Output(&outputA);
        edt.Run();		  
        edt.Input (& inputB);
        edt.Output(&outputB);
        edt.Run();

        for (t = 0 ; t < _tinput.GetT(); t++) {
          for (z = 0; z < _tinput.GetZ(); z++) {
            for (y = 0; y < _tinput.GetY(); y++) {
              for (x = 0; x < _tinput.GetX(); x++) {
                _dmap(x, y, z, t) = sqrt(outputA(x, y, z, t)) - sqrt(outputB(x, y, z, t));
              }
            }
          }
        }
      }

      // Linear Interpolate Dmap _dmap to _rdmap
      {
        double i,j,k,l;
        LinearInterpolateImageFunction interpolator;
        interpolator.Input(&_dmap);
        interpolator.Initialize();

        for (t = 0 ; t < _rdmap.GetT(); t++) {
          for (z = 0; z < _rdmap.GetZ(); z++) {
            for (y = 0; y < _rdmap.GetY(); y++) {
              for (x = 0; x < _rdmap.GetX(); x++) {
                i = x; j = y; k = z; l = t;
                _rdmap.ImageToWorld(i,j,k);
                _dmap.WorldToImage(i,j,k);
                _rdmap.PutAsDouble(x,y,z,t,interpolator.Evaluate(i,j,k,l));
              }
            }
          }
        }
      }

      // Put value back to Resampled Image _rinput < 0 if _rinput == 0 let the neighbor vote.
      {
        for (t = 0 ; t < _rdmap.GetT(); t++) {
          for (z = 0; z < _rdmap.GetZ(); z++) {
            for (y = 0; y < _rdmap.GetY(); y++) {
              for (x = 0; x < _rdmap.GetX(); x++) {
                dcurrent = _rdmap.GetAsDouble(x,y,z,t);
                if (dcurrent < 0 && current >= _rinput.GetAsDouble(x,y,z,t) && dcurrent < _rcdmap.GetAsDouble(x,y,z,t)) {
                  _rcdmap.PutAsDouble(x,y,z,t,dcurrent);
                } else if (dcurrent <= 0 && current >= _rinput.GetAsDouble(x,y,z,t) && dcurrent < _rcdmap.GetAsDouble(x,y,z,t)) {
                  sum = 0, sumcount = 0;
                  if (x > 0) {
                    sum += _rdmap.GetAsDouble(x-1,y,z,t); 
                    sumcount ++;
                  }
                  if (x < _rinput.GetX() - 1) {
                    sum += _rdmap.GetAsDouble(x+1,y,z,t); 
                    sumcount ++;
                  }
                  if (y > 0) {
                    sum += _rdmap.GetAsDouble(x,y-1,z,t); 
                    sumcount ++;
                  } 
                  if (y < _rinput.GetY() - 1) {
                    sum += _rdmap.GetAsDouble(x,y+1,z,t); 
                    sumcount++;
                  }
                  if (z > 0) {
                    sum += _rdmap.GetAsDouble(x,y,z-1,t); 
                    sumcount++;
                  }
                  if (z < _rinput.GetZ() - 1) {
                    sum += _rdmap.GetAsDouble(x,y,z+1,t); 
                    sumcount++;
                  }
                  sum = sum/sumcount;
                  if (sum <= 0) {
                    _rcdmap.PutAsDouble(x,y,z,t,dcurrent);
                  }
                } else if (dcurrent > 0 && _rinput.GetAsDouble(x,y,z,t) <= current && dcurrent < _rcdmap.GetAsDouble(x,y,z,t)) {
                  _rinput.PutAsDouble(x,y,z,t,current);
                  _rcdmap.PutAsDouble(x,y,z,t,dcurrent);
                } else if (dcurrent >= 0 && _rinput.GetAsDouble(x,y,z,t) <= current && dcurrent < _rcdmap.GetAsDouble(x,y,z,t)) {
                  if (sum > 0) {
                    _rinput.PutAsDouble(x,y,z,t,current);
                    _rcdmap.PutAsDouble(x,y,z,t,dcurrent);
                  }
                }
              }
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
  int    new_x, new_y, new_z, x, y, z, t,labelcount;
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
  new_x = int(Input()->GetX() * old_xsize / size);
  new_y = int(Input()->GetY() * old_ysize / size);
  new_z = int(Input()->GetZ() * old_zsize / size);

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
  for (current = min; current <= max; current++) {
    // Threshold
    sumcount = 0;
    for (t = 0; t < _tinput.GetT(); t++){
      for (z = 0; z < _tinput.GetZ(); z++) {
	      for (y = 0; y < _tinput.GetY(); y++) {
	        for (x = 0; x < _tinput.GetX(); x++) {
	          if (Input()->GetAsDouble(x, y, z, t) < current /* || Input()->GetAsDouble(x, y, z, t) > current */) {
	            _tinput(x, y, z, t) = 0;
            } else {
              _tinput(x, y, z, t) = 1;
            }
            if (Input()->GetAsDouble(x, y, z, t) >= current &&
                Input()->GetAsDouble(x, y, z, t) <  current + 1) {
              sumcount ++;
            }
          }
        }
      }
    }

    if (iround(current) % 20 == 0) {
      if (max < current + 19) {
        cout << "Doing outside DT for value >= : "<< current << " to " << max << endl;
        cout << "Doing inside DT for value >= : "<< current << " to " << max << endl;
      } else {
        cout << "Doing outside DT for value >= : "<< current << " to " << current+19 << endl;
        cout << "Doing inside DT for value >= : "<< current << " to " << current+19 << endl;
      }
    }

    if (sumcount > 0) {
      labelcount ++;
      // Dmap _tinput to _dmap
      {
        RealImage inputA, inputB, outputA, outputB;

        // Default mode
        DistanceTransformType edt(DistanceTransformType::DT_3D);

        // Threshold image
        inputA = _tinput;
        inputB = _tinput;
        for (t = 0; t < _tinput.GetT(); t++) {
          for (z = 0; z < _tinput.GetZ(); z++) {
            for (y = 0; y < _tinput.GetY(); y++) {
              for (x = 0; x < _tinput.GetX(); x++) {
                if (_tinput(x, y, z, t) > 0.5) {
                  inputA(x, y, z, t) = 1;
                  inputB(x, y, z, t) = 0;
                } else {
                  inputA(x, y, z, t) = 0;
                  inputB(x, y, z, t) = 1;
                }
              }
            }
          }
        }

        // Calculate EDT
        edt.Input (& inputA);
        edt.Output(&outputA);
        edt.Run();		  
        edt.Input (& inputB);
        edt.Output(&outputB);
        edt.Run();

        for (t = 0 ; t < _tinput.GetT(); t++) {
          for (z = 0; z < _tinput.GetZ(); z++) {
            for (y = 0; y < _tinput.GetY(); y++) {
              for (x = 0; x < _tinput.GetX(); x++) {
                _dmap(x, y, z, t)  = sqrt(outputA(x, y, z, t)) - sqrt(outputB(x, y, z, t));
              }
            }
          }
        }
      }

      // Linear Interpolate Dmap _dmap to _rdmap
      {
        double i,j,k,l;
        LinearInterpolateImageFunction interpolator;
        interpolator.Input(&_dmap);
        interpolator.Initialize();

        for (t = 0 ; t < _rdmap.T(); t++) {
          for (z = 0; z < _rdmap.Z(); z++) {
            for (y = 0; y < _rdmap.Y(); y++) {
              for (x = 0; x < _rdmap.X(); x++) {
                i = x; j = y; k = z; l = t;
                _rdmap.ImageToWorld(i,j,k);
                _dmap.WorldToImage(i,j,k);
                _rdmap.PutAsDouble(x,y,z,t,interpolator.Evaluate(i,j,k,l));
              }
            }
          }
        }
      }

      // Put value back to Resampled Image _rinput < 0 if _rinput == 0 let the neighbor vote.
      {
        for (t = 0 ; t < _rdmap.T(); t++) {
          for (z = 0; z < _rdmap.Z(); z++) {
            for (y = 0; y < _rdmap.Y(); y++) {
              for (x = 0; x < _rdmap.X(); x++) {
                if (_rdmap.GetAsDouble(x, y, z, t) < 0) {
                  _rinput.PutAsDouble(x, y, z, t, current);
                } else if (_rdmap.GetAsDouble(x, y, z, t) <= 0) {
                  sum = 0; sumcount = 0;
                  if (x > 0) {
                    sum += _rdmap.GetAsDouble(x-1,y,z,t); 
                    sumcount ++;
                  }
                  if (x < _rinput.GetX() - 1) {
                    sum += _rdmap.GetAsDouble(x+1,y,z,t); 
                    sumcount ++;
                  }
                  if (y > 0) {
                    sum += _rdmap.GetAsDouble(x,y-1,z,t); 
                    sumcount ++;
                  } 
                  if (y < _rinput.GetY() - 1) {
                    sum += _rdmap.GetAsDouble(x,y+1,z,t); 
                    sumcount++;
                  }
                  if (z > 0) {
                    sum += _rdmap.GetAsDouble(x,y,z-1,t); 
                    sumcount++;
                  }
                  if (z < _rinput.GetZ() - 1) {
                    sum += _rdmap.GetAsDouble(x,y,z+1,t); 
                    sumcount++;
                  }
                  sum = sum/sumcount;
                  if (sum <= 0) {
                    _rinput.PutAsDouble(x,y,z,t,current);
                  }
                }
              }
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
