/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
 * Copyright 2013-2017 Andreas Schuh
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

#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/IOConfig.h"
#include "mirtk/GenericImage.h"
#include "mirtk/VoxelFunction.h"
#include "mirtk/Transformations.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <target> <output> <dof> [options]" << endl;
  cout << "       " << name << " <target> <dof> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Computes the Jacobian determinant of the given transformation at" << endl;
  cout << "  each target image voxel. If no <output> image name is given," << endl;
  cout << "  the statistics of the Jacobian determinant distribution are printed" << endl;
  cout << "  to STDOUT instead, such as min, max, and mean." << endl;
  cout << endl;
  cout << "  By default, the output Jacobian determinant map has integral voxel type" << endl;
  cout << "  with Jacobian determinant values scaled by 100 prior to casting." << endl;
  cout << "  Use either option :option:`-float` or :option:`-double` to output the" << endl;
  cout << "  unscaled determinants." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  target   Target image." << endl;
  cout << "  output   Output Jacobian map." << endl;
  cout << "  dof      Input transformation." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -total               Total Jacobian. (default)" << endl;
  cout << "  -local               Local Jacobian only." << endl;
  cout << "  -global              Global Jacobian only." << endl;
  cout << "  -relative            Local Jacobian divided by global Jacobian." << endl;
  cout << "  -log                 Log of total Jacobian." << endl;
  cout << "  -padding <value>     Target padding value. (default: none)" << endl;
  cout << "  -outside <value>     Outside determinant value. (default: 0)" << endl;
  cout << "  -Tt <value>          Temporal origin of target image. (default: _torigin)" << endl;
  cout << "  -St <value>          Temporal origin of source image. (default: _torigin + _dt)" << endl;
  cout << "  -ss [fast]           Force use of scaling-and-squaring for SV FFD. (default: SV FFD integration method)" << endl;
  cout << "  -noss                Do not use scaling-and-squaring for SV FFD." << endl;
  cout << "  -disp [on|off]       Evaluate Jacobian of linearly interpolated displacement field. (default: off)" << endl;
  cout << "  -velo [on|off]       Evaluate Jacobian of velocity field parameterization for SV/TD FFD. (default: off)" << endl;
  cout << "  -float               Output values as single-precision floating point." << endl;
  cout << "  -double              Output values as double-precision floating point." << endl;
  cout << "  -real                Output values as floating point with default precision." << endl;
  cout << "  -threshold <value>   Lower Jacobian determinant threshold used when computing -log. (default: 1e-4)" << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

/// Enumeration of possible Jacobian determinant output values
enum JacobianMode {
  DefaultJacobian,
  LocalJacobian,
  GlobalJacobian,
  RelativeJacobian,
  TotalJacobian,
  LogJacobian,
  AbsLogJacobian,
  TotalAndLogJacobian
};

// -----------------------------------------------------------------------------
struct JacobianImpl : public VoxelReduction
{
  GreyImage            *_image;
  BaseImage            *_jacobian;
  BinaryImage          *_mask;
  int                   _nvox;
  const Transformation *_dof;
  double               _t, _t0;
  JacobianMode         _mode;
  int                  _padding;
  double               _outside;
  int                  _m, _n;
  int                  _noss;
  int                  _dtype;
  double               _threshold;

  JacobianImpl()
  :
    _image   (nullptr),
    _jacobian(nullptr),
    _mask    (nullptr),
    _nvox    (0),
    _dof     (nullptr),
    _t       (NaN),
    _t0      (NaN),
    _mode    (DefaultJacobian),
    _padding (MIN_GREY),
    _outside (.0),
    _m       (0),
    _n       (0),
    _noss    (-1),
    _dtype   (MIRTK_VOXEL_GREY),
    _threshold(.0001)
  {}

  JacobianImpl(const JacobianImpl &other)
  :
    _image   (other._image),
    _jacobian(other._jacobian),
    _mask    (other._mask),
    _nvox    (other._nvox),
    _dof     (other._dof),
    _t       (other._t),
    _t0      (other._t0),
    _mode    (other._mode),
    _padding (other._padding),
    _outside (other._outside),
    _m       (other._m),
    _n       (other._n),
    _noss    (other._noss),
    _dtype   (other._dtype),
    _threshold(other._threshold)
  {}

  void split(const JacobianImpl &other)
  {
    _m = _n = 0;
  }

  void join(const JacobianImpl &other)
  {
    _m += other._m;
    _n += other._n;
  }

  double Jacobian(int i, int j, int k, const BinaryPixel *mask)
  {
    if (*mask == 0) {
      return _outside;
    }

    ++_m;

    double jac = numeric_limits<double>::quiet_NaN();

    double x = i, y = j, z = k;
    _image->ImageToWorld(x, y, z);
    switch (_mode) {
      case DefaultJacobian:
      case TotalJacobian: {
        jac = _dof->Jacobian(x, y, z, _t, _t0);
        if (jac < .0) ++_n;
      } break;
      case LocalJacobian: {
        jac = _dof->LocalJacobian(x, y, z, _t, _t0);
        if (jac < .0) ++_n;
      }break;
      case GlobalJacobian: {
        jac = _dof->GlobalJacobian(x, y, z, _t, _t0);
      }break;
      case RelativeJacobian: {
        jac = _dof->LocalJacobian(x, y, z, _t, _t0);
        if (jac < .0) ++_n;
        jac /= _dof->GlobalJacobian(x, y, z, _t, _t0);
      } break;
      case LogJacobian: {
        jac = _dof->Jacobian(x, y, z, _t, _t0);
        if (jac < .0) ++_n;
        jac = log(max(jac, _threshold));
      } break;
      case AbsLogJacobian: {
        jac = _dof->Jacobian(x, y, z, _t, _t0);
        if (jac < .0) ++_n;
        jac = fabs(log(max(jac, _threshold)));
      } break;
      case TotalAndLogJacobian: {
        jac = _dof->Jacobian(x, y, z, _t, _t0);
        if (jac < .0) ++_n;
      } break;
    }

    return IsNaN(jac) ? _outside : jac;
  }

  void operator()(int i, int j, int k, int, const BinaryPixel *mask, GreyPixel *out)
  {
    double jac = Jacobian(i, j, k, mask);
    *out = static_cast<GreyPixel>(100.0 * jac);
    if (_mode == TotalAndLogJacobian) {
      *(out + _nvox) = static_cast<GreyPixel>(100.0 * log(max(jac, _threshold)));
    }
  }

  void operator()(int i, int j, int k, int, const BinaryPixel *mask, float *out)
  {
    double jac = Jacobian(i, j, k, mask);
    *out = static_cast<float>(100.0 * jac);
    if (_mode == TotalAndLogJacobian) {
      *(out + _nvox) = static_cast<float>(100.0 * log(max(jac, _threshold)));
    }
  }

  void operator()(int i, int j, int k, int, const BinaryPixel *mask, double *out)
  {
    double jac = Jacobian(i, j, k, mask);
    *out = 100.0 * jac;
    if (_mode == TotalAndLogJacobian) {
      *(out + _nvox) = 100.0 * log(max(jac, _threshold));
    }
  }

  void Run()
  {
    _nvox = _image->NumberOfSpatialVoxels();
    ImageAttributes attr = _image->Attributes();
    attr._t = 1, attr._dt = .0;
    // Default arguments
    if (IsNaN(_t0)) _t0 = _image->ImageToTime(0);
    if (IsNaN(_t )) _t  = _t0 + _image->GetTSize();
    // Initialize output image
    if (_mode == TotalAndLogJacobian) _jacobian->Initialize(attr, 2);
    else                              _jacobian->Initialize(attr, 1);
    // Use scaling and squaring in case of SV FFD
    const BSplineFreeFormTransformationSV *svffd = nullptr;
    if (_mode != GlobalJacobian && _noss != 0) {
      svffd = dynamic_cast<const BSplineFreeFormTransformationSV *>(_dof);
      if (!svffd) {
        const MultiLevelTransformation *mffd;
        mffd = dynamic_cast<const MultiLevelTransformation *>(_dof);
        if (mffd && mffd->NumberOfLevels() == 1 && mffd->GetGlobalTransformation()->IsIdentity()) {
          svffd = dynamic_cast<const BSplineFreeFormTransformationSV *>(mffd->GetLocalTransformation(0));
        }
      }
    }
    if (svffd && (_noss > 0 || svffd->IntegrationMethod() == FFDIM_SS || svffd->IntegrationMethod() == FFDIM_FastSS)) {
      const FFDIntegrationMethod ffdim = svffd->IntegrationMethod();
      if (_noss > 0) {
        if (_noss == 2) {
          const_cast<BSplineFreeFormTransformationSV *>(svffd)->IntegrationMethod(FFDIM_FastSS);
        } else {
          const_cast<BSplineFreeFormTransformationSV *>(svffd)->IntegrationMethod(FFDIM_SS);
        }
      }
#if 0
      // Compute all partial derivatives using scaling and squaring
      if (verbose) cout << "Computing Jac using scaling and squaring" << endl;
      GenericImage<double> d;
      MIRTK_START_TIMING();
      svffd->ScalingAndSquaring<double>(attr, nullptr, &d, nullptr, nullptr, svffd->UpperIntegrationLimit(_t, _t0));
      MIRTK_DEBUG_TIMING(1, "computation of Jac using scaling and squaring");
      // Compute determinants of Jacobian matrices
      if (verbose) cout << "Computing det(Jac)" << endl;
      MIRTK_RESET_TIMING();
      const int xx = 0 * _nvox, xy = 1 * _nvox, xz = 2 * _nvox;
      const int yx = 3 * _nvox, yy = 4 * _nvox, yz = 5 * _nvox;
      const int zx = 6 * _nvox, zy = 7 * _nvox, zz = 8 * _nvox;
      double jac;
      for (int idx = 0; idx < _nvox; ++idx) {
        if (_mask->Get(idx) != 0) {
          jac = d(xx + idx) * d(yy + idx) * d(zz + idx)
              + d(xy + idx) * d(yz + idx) * d(zx + idx)
              + d(xz + idx) * d(yx + idx) * d(zy + idx)
              - d(xz + idx) * d(yy + idx) * d(zx + idx)
              - d(xy + idx) * d(yx + idx) * d(zz + idx)
              - d(xx + idx) * d(yz + idx) * d(zy + idx);
          _m++; if (jac < .0) _n++;
          if (_mode == LogJacobian) {
            jac = fabs(log(max(jac, _threshold)));
          } else if (_mode == TotalAndLogJacobian) {
            _jacobian->PutAsDouble(_nvox + idx, 100.0 * log(max(jac, _threshold)));
          }
          jac *= 100.;
        } else {
          jac = _outside;
        }
        _jacobian->PutAsDouble(idx, jac);
      }
      MIRTK_DEBUG_TIMING(1, "computation of det(Jac)");
#else
      // Compute Jacobian using scaling and squaring
      const double T = svffd->UpperIntegrationLimit(_t, _t0);
      switch (_mode) {
        case LogJacobian: {
          if (verbose) cout << "Computing log(det(Jac)) using scaling and squaring" << endl;
          GenericImage<double> lj;
          MIRTK_START_TIMING();
          svffd->ScalingAndSquaring<double>(attr, nullptr, nullptr, nullptr, &lj, T);
          MIRTK_DEBUG_TIMING(1, "computation of log(det(Jac)) using scaling and squaring");
          for (int idx = 0; idx < _nvox; ++idx) {
            if (_mask->Get(idx) != 0) {
              _jacobian->PutAsDouble(idx, 100.0 * lj(idx));
              ++_m;
            } else {
              _jacobian->PutAsDouble(idx, _outside);
            }
          }
        } break;
        case AbsLogJacobian: {
          if (verbose) cout << "Computing log(det(Jac)) using scaling and squaring" << endl;
          GenericImage<double> lj;
          MIRTK_START_TIMING();
          svffd->ScalingAndSquaring<double>(attr, nullptr, nullptr, nullptr, &lj, T);
          MIRTK_DEBUG_TIMING(1, "computation of log(det(Jac)) using scaling and squaring");
          for (int idx = 0; idx < _nvox; ++idx) {
            if (_mask->Get(idx) != 0) {
              _jacobian->PutAsDouble(idx, 100.0 * fabs(lj(idx)));
              ++_m;
            } else {
              _jacobian->PutAsDouble(idx, _outside);
            }
          }
        } break;
        case TotalAndLogJacobian: {
          if (verbose) cout << "Computing det(Jac) and its log using scaling and squaring" << endl;
          GenericImage<double> dj, lj;
          MIRTK_START_TIMING();
          svffd->ScalingAndSquaring<double>(attr, nullptr, nullptr, &dj, &lj, T);
          MIRTK_DEBUG_TIMING(1, "computation of det(Jac) and its log using scaling and squaring");
          for (int idx = 0; idx < _nvox; ++idx) {
            if (_mask->Get(idx) != 0) {
              _jacobian->PutAsDouble(idx, 100.0 * dj(idx));
              _jacobian->PutAsDouble(idx + _nvox, 100.0 * lj(idx));
              ++_m, _n += (dj(idx) < .0 ? 1 : 0);
            } else {
              _jacobian->PutAsDouble(idx, _outside);
            }
          }
        } break;
        default: {
          if (verbose) cout << "Computing det(Jac) using scaling and squaring" << endl;
          GenericImage<double> dj;
          MIRTK_START_TIMING();
          svffd->ScalingAndSquaring<double>(attr, nullptr, nullptr, &dj, nullptr, T);
          MIRTK_DEBUG_TIMING(1, "computation of det(Jac) using scaling and squaring");
          for (int idx = 0; idx < _nvox; ++idx) {
            if (_mask->Get(idx) != 0) {
              _jacobian->PutAsDouble(idx, 100.0 * dj(idx));
              ++_m, _n += (dj(idx) < .0 ? 1 : 0);
            } else {
              _jacobian->PutAsDouble(idx, _outside);
            }
          }
        } break;
      }
#endif
      const_cast<BSplineFreeFormTransformationSV *>(svffd)->IntegrationMethod(ffdim);
    // Otherwise, evaluate Jacobian for each voxel separately
    } else {
      if (verbose) cout << "Computing Jacobian of non-SV FFD (or -noss option given)" << endl;
      GreyImage            *iout = nullptr;
      GenericImage<float>  *fout = nullptr;
      GenericImage<double> *dout = nullptr;
      MIRTK_START_TIMING();
      if ((iout = dynamic_cast<GreyImage *>(_jacobian))) {
        ParallelForEachVoxel(attr, _mask, iout, *this);
      } else if ((fout = dynamic_cast<GenericImage<float> *>(_jacobian))) {
        ParallelForEachVoxel(attr, _mask, fout, *this);
      } else if ((dout = dynamic_cast<GenericImage<double> *>(_jacobian))) {
        ParallelForEachVoxel(attr, _mask, dout, *this);
      } else {
        cerr << "Output image data type must be either GreyPixel, float, or double" << endl;
        exit(1);
      }
      MIRTK_DEBUG_TIMING(1, "computation of Jacobian determinant");
    }
  }
};

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  REQUIRES_POSARGS(2);

  const char *target_name = POSARG(1);
  const char *output_name = nullptr;
  const char *dofin_name  = nullptr;
  if      (NUM_POSARGS == 2) dofin_name = POSARG(2), verbose     = 1;
  else if (NUM_POSARGS == 3) dofin_name = POSARG(3), output_name = POSARG(2);
  else {
    FatalError("Invalid number of positional arguments!");
  }

  InitializeIOLibrary();
  UniquePtr<GreyImage> target(new GreyImage(target_name));
  UniquePtr<Transformation> dof(Transformation::New(dofin_name));

  JacobianImpl impl;
  bool asdisp = false;
  bool asvelo = false;

  for (ALL_OPTIONS) {
    if (OPTION("-padding")) {
      PARSE_ARGUMENT(impl._padding);
    }
    else if (OPTION("-outside")) {
      PARSE_ARGUMENT(impl._outside);
    }
    else if (OPTION("-threshold")) {
      PARSE_ARGUMENT(impl._threshold);
    }
    else if (OPTION("-local"))     impl._mode = LocalJacobian;
    else if (OPTION("-global"))    impl._mode = GlobalJacobian;
    else if (OPTION("-total"))     impl._mode = TotalJacobian;
    else if (OPTION("-log"))       impl._mode = LogJacobian;
    else if (OPTION("-abslog"))    impl._mode = AbsLogJacobian;
    else if (OPTION("-total+log")) impl._mode = TotalAndLogJacobian;
    else if (OPTION("-relative"))  impl._mode = RelativeJacobian;
    else if (OPTION("-fluid")) ; // ignore deprecated/obsolete option
    else if (OPTION("-Tt")) {
      PARSE_ARGUMENT(impl._t0);
    }
    else if (OPTION("-St")) {
      PARSE_ARGUMENT(impl._t);
    }
    else if (OPTION("-ss")) {
      impl._noss = 1;
      if (HAS_ARGUMENT) {
        bool ss;
        const string arg = ARGUMENT;
        if (ToLower(arg) == "fast") {
          impl._noss = 2;
        } else if (FromString(arg, ss)) {
          impl._noss = (ss ? 1 : 0);
        } else {
          FatalError("Invalid -ss option argument: " << arg);
        }
      }
    }
    else if (OPTION("-disp")) {
      asdisp = true;
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(asdisp);
      }
    }
    else if (OPTION("-velo")) {
      asvelo = true;
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(asvelo);
      }
    }
    else if (OPTION("-noss"))   impl._noss  = 0;
    else if (OPTION("-float"))  impl._dtype = MIRTK_VOXEL_FLOAT;
    else if (OPTION("-double")) impl._dtype = MIRTK_VOXEL_DOUBLE;
    else if (OPTION("-real"))   impl._dtype = MIRTK_VOXEL_REAL;
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (asvelo) {
    BSplineFreeFormTransformationSV *svffd = nullptr;
    BSplineFreeFormTransformationTD *tdffd = nullptr;
    svffd = dynamic_cast<BSplineFreeFormTransformationSV *>(dof.get());
    tdffd = dynamic_cast<BSplineFreeFormTransformationTD *>(dof.get());
    if (svffd) {
      svffd->ScaleVelocities(svffd->T() / svffd->NumberOfSteps());
      dof.reset(new BSplineFreeFormTransformation3D(*svffd));
    } else if (tdffd) {
      dof.reset(new BSplineFreeFormTransformation4D(*tdffd));
    } else {
      FatalError("Input must be SV/TD FFD when -velo option is used");
    }
  }
  if (asdisp && impl._mode != GlobalJacobian) {
    if (impl._mode == LocalJacobian) {
      FatalError("Option -disp not implemented for -local Jacobian");
    }
    GenericImage<double> disp(target->Attributes(), 3);
    dof->Displacement(disp, impl._t, impl._t0);
    dof.reset(new LinearFreeFormTransformation3D(disp));
  }

  impl._image    = target.get();
  impl._jacobian = target.get();
  impl._dof      = dof.get();

  const int nvox = target->NumberOfSpatialVoxels();

  BinaryImage mask(target->Attributes(), 1);
  for (int idx = 0; idx < nvox; ++idx) {
    mask(idx) = (target->GetAsDouble(idx) > impl._padding ? 1 : 0);
  }
  impl._mask = &mask;

  UniquePtr<BaseImage> output;
  if (impl._dtype != target->GetDataType()) {
    output.reset(BaseImage::New(impl._dtype));
    impl._jacobian = output.get();
  }
  impl.Run();
  if (impl._jacobian == impl._image) {
    output.reset(target.release());
  }

  if (verbose) {
    int    num = 0;
    double jac, min = inf, max = -inf, avg = 0.;
    for (int n = 0; n < nvox; ++n) {
      if (mask(n) != 0) {
        jac = output->GetAsDouble(n);
        if (jac < min) min = jac;
        if (jac > max) max = jac;
        avg += jac;
        num += 1;
      }
    }
    min /= 100.0, max /= 100.0, avg /= 100.0;
    if (num > 0) avg /= num;
    cout << "Minimum value = " << min << endl;
    cout << "Maximum value = " << max << endl;
    cout << "Average value = " << avg << endl;
    if (impl._n > 0) {
      cout << endl;
      cout << "Number of voxels with negative Jacobian determinant = ";
      cout << impl._n << " (" << 100.0 * impl._n / impl._m << "%)" << endl;
    }
  }

  if (impl._dtype == MIRTK_VOXEL_FLOAT) {
    GenericImage<float> *jacobian;
    jacobian = dynamic_cast<GenericImage<float> *>(output.get());
    *jacobian /= 100.0f;
  } else if (impl._dtype == MIRTK_VOXEL_DOUBLE) {
    GenericImage<double> *jacobian;
    jacobian = dynamic_cast<GenericImage<double> *>(output.get());
    *jacobian /= 100.0;
  }

  if (output_name) {
    output->Write(output_name);
  }
  return 0;
}
