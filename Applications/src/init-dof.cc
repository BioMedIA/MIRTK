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

#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/IOConfig.h"
#include "mirtk/PointSamples.h"

#include "mirtk/Transformation.h"
#include "mirtk/HomogeneousTransformation.h"
#include "mirtk/FreeFormTransformation.h"
#include "mirtk/FreeFormTransformation3D.h"
#include "mirtk/MultiLevelTransformation.h"

#include "mirtk/RigidTransformation.h"
#include "mirtk/SimilarityTransformation.h"
#include "mirtk/AffineTransformation.h"
#include "mirtk/BSplineFreeFormTransformation3D.h"
#include "mirtk/BSplineFreeFormTransformationSV.h"
#include "mirtk/MultiLevelFreeFormTransformation.h"

#if defined(HAVE_VTK) && MIRTK_IO_WITH_VTK
#  include "vtkSmartPointer.h"
#  include "vtkPointSet.h"
#  include "vtkPointData.h"
#  include "vtkDataArray.h"
#  include "mirtk/PointSetIO.h"
#endif

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <dofout> [output option] [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  This tool either creates a new affine transformation with the given parameters" << endl;
  cout << "  or approximates an affine transformation or non-rigid deformation given" << endl;
  cout << "  a set of corresponding landmarks or point displacements (see :option:`-displacements`)." << endl;
  cout << endl;
  cout << "  The output transformation is by default an affine transformation." << endl;
  cout << "  One of the output options below can be used to request a different" << endl;
  cout << "  type of output transformation. This option should be given right" << endl;
  cout << "  after the output transformation name as it may be altered using" << endl;
  cout << "  one of the other output related options (e.g. disable :option:`-shearing`)." << endl;
  cout << endl;
  cout << "  An affine transformation is the composition of a shearing followed by" << endl;
  cout << "  a scaling, a rotation, and a translation. The homogeneous transformation" << endl;
  cout << "  matrix, A, is given by the matrix product:" << endl;
  cout << endl;
  cout << "  .. math::" << endl;
  cout << endl;
  cout << "     A = Translate * Rotate * Scale * Shear, where Rotate = (Rz Ry Rx)^T" << endl;
  cout << endl;
  cout << "Output options:" << endl;
  cout << "  -rigid         Output rigid transformation." << endl;
  cout << "  -similarity    Output similarity transformation." << endl;
  cout << "  -affine        Output affine transformation." << endl;
  cout << "  -affine-mffd   Output affine transformation as global component of MFFD." << endl;
  cout << "  -mffd          Output multi-level free-form deformation (MFFD)." << endl;
  cout << "  -ffd           Output free-form deformation (FFD)." << endl;
  cout << "  -svffd         Output diffeomorphic free-form deformation (SV FFD)." << endl;
  cout << endl;
  cout << "  -[no]translations   Allow/disallow translations." << endl;
  cout << "  -[no]rotations      Allow/disallow rotations." << endl;
  cout << "  -[no]scaling        Allow/disallow scaling." << endl;
  cout << "  -[no]shearing       Allow/disallow shearing." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -dofin <dof>      Read affine parameters from input transformation." << endl;
  cout << endl;
  cout << "  -orientation <x1> <x2> <x3>  <y1> <y2> <y3>  <z1> <z2> <z3>" << endl;
  cout << "      Set upper 3x3 (rotation) matrix entries" << endl;
  cout << "      (e.g., image orientation from image info output)." << endl;
  cout << endl;
  cout << "  -image-reset <file>" << endl;
  cout << "      Set transformation to difference between spatial image to world mapping" << endl;
  cout << "      of given image and the default lattice with axes parallel to the axes of" << endl;
  cout << "      the world coordinate system and image center at the origin. The affine" << endl;
  cout << "      output transformation maps points of the unoriented world system with" << endl;
  cout << "      origin at the image center to the world coordinates of the input image." << endl;
  cout << "  -image-rotation <file>" << endl;
  cout << "      Set transformation to difference between spatial image to world mapping" << endl;
  cout << "      of given image and an image lattice with axes parallel to the axes of" << endl;
  cout << "      the world coordinate system. The affine output transformation maps points" << endl;
  cout << "      of the unoriented world system to the world coordinates of the input image." << endl;
  cout << "      Hence, the domain of the map is the same as the world domain of the image," << endl;
  cout << "      except that the orientation matrix is equal the identity matrix." << endl;
  cout << endl;
  cout << "  -tx <tx>          Translation along x axis in mm.             (default: 0)" << endl;
  cout << "  -ty <ty>          Translation along y axis in mm.             (default: 0)" << endl;
  cout << "  -tz <tz>          Translation along z axis in mm.             (default: 0)" << endl;
  cout << "  -rx <rx>          Rotation around x axis in degrees.          (default: 0)" << endl;
  cout << "  -ry <ry>          Rotation around x axis in degrees.          (default: 0)" << endl;
  cout << "  -rz <rz>          Rotation around x axis in degrees.          (default: 0)" << endl;
  cout << "  -sx <sx>          Scaling of x axis in percentage.            (default: 100)" << endl;
  cout << "  -sy <sy>          Scaling of y axis in percentage.            (default: 100)" << endl;
  cout << "  -sz <sz>          Scaling of z axis in percentage.            (default: 100)" << endl;
  cout << "  -s <s>            Isotropic scaling in percentage.            (default: 100)" << endl;
  cout << "  -sxy <sxy>        Skew angle between x and y axis in degrees. (default: 0)" << endl;
  cout << "  -sxz <sxz>        Skew angle between x and z axis in degrees. (default: 0)" << endl;
  cout << "  -syz <syz>        Skew angle between y and z axis in degrees. (default: 0)" << endl;
  cout << endl;
  cout << "  -target <image>   Target image. Required when output is a FFD. (default: none)" << endl;
  cout << "  -source <image>   Source image. Used to initialize translation. (default: none)" << endl;
  cout << "  -Tp <value>       Target background value. (default: none)" << endl;
  cout << "  -Sp <value>       Source background value. (default: none)" << endl;
  cout << endl;
  cout << "  -align-top        Align top-most slices of :option:`-target` and :option:`-source` image. (default: center)" << endl;
  cout << endl;
  cout << "  -displacements, -disp <pset1> <pset2> [<cor12>]" << endl;
  cout << "      Create transformation which minimizes the mean squared distance between pairs" << endl;
  cout << "      of fiducial point sets <pset1> and <pset2>. When a single <pset> is given," << endl;
  cout << "      it creates a transformation which approximates the displacement vectors stored" << endl;
  cout << "      in the \"vectors\" point set attributes array. Option can be used multiple times." << endl;
  cout << endl;
  cout << "  -approximate <dofin>              Create transformation which approximates" << endl;
  cout << "                                    another given transformation." << endl;
  cout << "                                    Option can be used multiple times." << endl;
  cout << endl;
  cout << "Free-form deformation options:" << endl;
  cout << "  -dx <dx>         Control point spacing of FFD in x. (default: target voxel size)" << endl;
  cout << "  -dy <dy>         Control point spacing of FFD in y. (default: target voxel size)" << endl;
  cout << "  -dz <dz>         Control point spacing of FFD in z. (default: target voxel size)" << endl;
  cout << "  -ds <ds>         Control point spacing of FFD.      (default: target voxel size)" << endl;
  cout << "  -subdiv <n>      Number of subdivisions to use for approximating a non-rigid deformation" << endl;
  cout << "                   from input :option:`-displacements`. (default: 4)" << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
/// Read correspondences from text file
Array<int> ReadCorrespondences(const char *fname, int n1, int n2)
{
  ifstream ifs(fname);
  if (!ifs.is_open()) {
    cerr << "Error: Failed to open correspondence input file: " << fname << endl;
    exit(1);
  }
  string line;
  int    t, s;
  Array<int> corr(n1, -1);
  while (getline(ifs, line)) {
    istringstream is(line);
    if (!(is >> t >> s)) {
      cerr << "Error: Failed to read correspondence map from file: " << fname << endl;
      exit(1);
    }
    if (t < 0 || t >= n1) {
      cerr << "Error: Invalid target index (" << t << ") in correspondence map: " << fname << endl;
      exit(1);
    }
    if (s < 0 || s >= n2) {
      cerr << "Error: Invalid source index (" << s << ") in correspondence map: " << fname << endl;
      exit(1);
    }
    corr[t] = s;
  }
  for (t = 0; t < n1; ++t) {
    if (corr[t] == -1) {
      cerr << "Error: Missing correspondence for target index " << t << endl;
      exit(1);
    }
  }
  ifs.close();
  return corr;
}

// -----------------------------------------------------------------------------
FreeFormTransformation3D *
BSplineFFD(const AffineTransformation &dof, const ImageAttributes &attr, double sx, double sy, double sz)
{
  BSplineFreeFormTransformation3D *ffd = new BSplineFreeFormTransformation3D(attr, sx, sy, sz);

  const int N = ffd->NumberOfCPs();
  double *dx = new double[N];
  double *dy = new double[N];
  double *dz = new double[N];

  int n = 0;
  for (int k = 0; k < ffd->GetZ(); ++k)
  for (int j = 0; j < ffd->GetY(); ++j)
  for (int i = 0; i < ffd->GetX(); ++i, ++n) {
    dx[n] = i, dy[n] = j, dz[n] = k;
    ffd->LatticeToWorld(dx[n], dy[n], dz[n]);
    dof .Displacement  (dx[n], dy[n], dz[n]);
  }

  ffd->Interpolate(dx, dy, dz);

  delete[] dx;
  delete[] dy;
  delete[] dz;

  return ffd;
}

// -----------------------------------------------------------------------------
FreeFormTransformation3D *
BSplineSVFFD(const AffineTransformation &dof, ImageAttributes attr, double sx, double sy, double sz)
{
  BSplineFreeFormTransformationSV *ffd = new BSplineFreeFormTransformationSV();

  ffd->NumberOfSteps(64);
  ffd->MaxScaledVelocity(.0);
  ffd->Initialize(attr, sx, sy, sz, &dof);

  if (verbose) {
    double error = .0, x, y, z;
    GenericImage<double> disp(attr, 3);
    ffd->Displacement(disp);
    for (int k = 0; k < disp.GetZ(); ++k)
    for (int j = 0; j < disp.GetY(); ++j)
    for (int i = 0; i < disp.GetX(); ++i) {
      x = i, y = j, z = k;
      disp.ImageToWorld(x, y, z);
      dof .Displacement(x, y, z);
      x -= disp(i, j, k, 0);
      y -= disp(i, j, k, 1);
      z -= disp(i, j, k, 2);
      error += x*x + y*y + z*z;
    }
    error = sqrt(error / disp.GetX() * disp.GetY() * disp.GetZ());
    cout << "RMS error of B-spline SV FFD approximation = " << error << "\n" << endl;
  }

  return ffd;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  InitializeIOLibrary();
 
  EXPECTS_POSARGS(1);

  // Number of point samples used for approximation of rigid/affine input transformation
  const int number_of_point_samples = 20;

  // Instantiate affine transformation
  AffineTransformation dof;

  struct PointSetPair
  {
    const char *target_name; // Target fiducial points
    const char *source_name; // Source fiducial points
    const char *corrin_name; // Correspondences between target and source points
    PointSetPair() : target_name(NULL), source_name(NULL), corrin_name(NULL) {}
  };

  // Parse command-line arguments
  const char *dofout_name    = POSARG(1);
  const char *target_name    = NULL;  // Target image/FFD lattice read from NIfTI header
  const char *source_name    = NULL;  // Source image
  double      target_padding = numeric_limits<double>::quiet_NaN();
  double      source_padding = numeric_limits<double>::quiet_NaN();
  bool        align_top      = false; // Align top slice of target and source
  bool        align_centers  = false; // Align image centers of target and source
  bool        align_cogs     = false; // Align foreground centers of mass of target and source
  double      sx             = .0;    // Control point spacing in x = voxel size by default
  double      sy             = .0;    // Control point spacing in y = voxel size by default
  double      sz             = .0;    // Control point spacing in z = voxel size by default
  bool        similarity     = false; // Isotropic scaling
  bool        multilevel     = false; // Create multi-level FFD
  bool        nonrigid       = false; // Create non-rigid transformation
  bool        diffeomorphic  = false; // B-spline FFD or affine transformation by default
  int         nsubdiv        = 4;     // No. of subdivisions to approximate non-linear deformation
  Array<PointSetPair> psets;

  for (ALL_OPTIONS) {
    // Output options
    if (OPTION("-rigid")) {
      multilevel = false, nonrigid = false, diffeomorphic = true;
      dof.AllowTranslations(true);
      dof.AllowRotations(true);
      dof.AllowScaling(false);
      dof.AllowShearing(false);
    }
    else if (OPTION("-similarity")) {
      similarity = true;
      dof.AllowTranslations(true);
      dof.AllowRotations(true);
      dof.AllowScaling(true);
      dof.AllowShearing(false);
    }
    else if (OPTION("-affine") || OPTION("-affine-mffd")) {
      multilevel    = (strstr(OPTNAME, "mffd") != NULL);
      nonrigid      = false;
      diffeomorphic = true;
      dof.AllowTranslations(true);
      dof.AllowRotations(true);
      dof.AllowScaling(true);
      dof.AllowShearing(true);
    }
    else if (OPTION("-mffd") || OPTION("-ffd") || OPTION("-svffd")) {
      multilevel    = (strstr(OPTNAME, "mffd") != NULL);
      nonrigid      = true;
      diffeomorphic = (strstr(OPTNAME, "sv") != NULL);
    }
    // Enable/disable affine transformation parameters
    // Note: Considered only for approximation of input displacements!
    else if (OPTION("-translations"  )) dof.AllowTranslations(true);
    else if (OPTION("-notranslations")) dof.AllowTranslations(false);
    else if (OPTION("-rotations"     )) dof.AllowRotations(true);
    else if (OPTION("-norotations"   )) dof.AllowRotations(false);
    else if (OPTION("-scaling"       )) dof.AllowScaling(true);
    else if (OPTION("-noscaling"     )) dof.AllowScaling(false);
    else if (OPTION("-shearing"      )) dof.AllowShearing(true);
    else if (OPTION("-noshearing"    )) dof.AllowShearing(false);
    // Initialize upper 3x3 (rotation) matrix
    else if (OPTION("-orientation")) {
      Matrix m(4, 4);
      m(3, 3) = 1.0;
      const char *arg;
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
          arg = ARGUMENT;
          if (arg[0] == ',' && arg[1] == '\0') arg = ARGUMENT;
          m(r, c) = atof(arg);
        }
      }
      dof.PutMatrix(m);
    }
    else if (OPTION("-image-reset")) {
      BinaryImage image(ARGUMENT);
      ImageAttributes attr;
      attr._x  = image.X();
      attr._y  = image.Y();
      attr._z  = image.Z();
      attr._t  = image.T();
      attr._dx = image.XSize();
      attr._dy = image.YSize();
      attr._dz = image.ZSize();
      attr._dt = image.TSize();
      attr._torigin = image.GetTOrigin();
      dof.PutMatrix(image.GetImageToWorldMatrix() * attr.GetWorldToImageMatrix());
    }
    else if (OPTION("-image-rotation")) {
      BinaryImage image(ARGUMENT);
      ImageAttributes attr = image.Attributes();
      attr._xaxis[0] = 1.;
      attr._xaxis[1] = 0.;
      attr._xaxis[2] = 0.;
      attr._yaxis[0] = 0.;
      attr._yaxis[1] = 1.;
      attr._yaxis[2] = 0.;
      attr._zaxis[0] = 0.;
      attr._zaxis[1] = 0.;
      attr._zaxis[2] = 1.;
      dof.PutMatrix(image.GetImageToWorldMatrix() * attr.GetWorldToImageMatrix());
    }
    // Read affine parameters from input transformation
    else if (OPTION("-dofin")) {
      Transformation *dofin = Transformation::New(ARGUMENT);
      RigidTransformation      *dof6  = NULL;
      SimilarityTransformation *dof7  = NULL;
      AffineTransformation     *dof12 = NULL;
      MultiLevelTransformation *mffd  = NULL;
      // Attention: Order of casts or if clauses below matters!
      (dof12 = dynamic_cast<AffineTransformation     *>(dofin)) ||
      (dof7  = dynamic_cast<SimilarityTransformation *>(dofin)) ||
      (dof6  = dynamic_cast<RigidTransformation      *>(dofin)) ||
      (mffd  = dynamic_cast<MultiLevelTransformation *>(dofin));
      if (mffd) dof12 = mffd->GetGlobalTransformation();
      if (dof12) {
        dof.PutTranslationX(dof12->GetTranslationX());
        dof.PutTranslationY(dof12->GetTranslationY());
        dof.PutTranslationZ(dof12->GetTranslationZ());
        dof.PutRotationX   (dof12->GetRotationX());
        dof.PutRotationY   (dof12->GetRotationY());
        dof.PutRotationZ   (dof12->GetRotationZ());
        dof.PutScaleX      (dof12->GetScaleX());
        dof.PutScaleY      (dof12->GetScaleY());
        dof.PutScaleZ      (dof12->GetScaleZ());
        dof.PutShearXY     (dof12->GetShearXY());
        dof.PutShearXZ     (dof12->GetShearXZ());
        dof.PutShearYZ     (dof12->GetShearYZ());
      } else if (dof7) {
        dof.PutTranslationX(dof7->GetTranslationX());
        dof.PutTranslationY(dof7->GetTranslationY());
        dof.PutTranslationZ(dof7->GetTranslationZ());
        dof.PutRotationX   (dof7->GetRotationX());
        dof.PutRotationY   (dof7->GetRotationY());
        dof.PutRotationZ   (dof7->GetRotationZ());
        dof.PutScale       (dof7->GetScale());
      } else if (dof6) {
        dof.PutTranslationX(dof6->GetTranslationX());
        dof.PutTranslationY(dof6->GetTranslationY());
        dof.PutTranslationZ(dof6->GetTranslationZ());
        dof.PutRotationX   (dof6->GetRotationX());
        dof.PutRotationY   (dof6->GetRotationY());
        dof.PutRotationZ   (dof6->GetRotationZ());
      } else {
        delete dofin;
        cerr << "Only a rigid, similarity, affine, or MFFD transformation can be used as -dofin" << endl;
        exit(1);
      }
      delete dofin;
    }
    // Affine transformation parameters
    else if (OPTION("-tx"))  dof.PutTranslationX(atof(ARGUMENT));
    else if (OPTION("-ty"))  dof.PutTranslationY(atof(ARGUMENT));
    else if (OPTION("-tz"))  dof.PutTranslationZ(atof(ARGUMENT));
    else if (OPTION("-rx"))  dof.PutRotationX   (atof(ARGUMENT));
    else if (OPTION("-ry"))  dof.PutRotationY   (atof(ARGUMENT));
    else if (OPTION("-rz"))  dof.PutRotationZ   (atof(ARGUMENT));
    else if (OPTION("-sx"))  dof.PutScaleX      (atof(ARGUMENT));
    else if (OPTION("-sy"))  dof.PutScaleY      (atof(ARGUMENT));
    else if (OPTION("-sz"))  dof.PutScaleZ      (atof(ARGUMENT));
    else if (OPTION("-s"))   dof.PutScale       (atof(ARGUMENT));
    else if (OPTION("-sxy")) dof.PutShearXY     (atof(ARGUMENT));
    else if (OPTION("-sxz")) dof.PutShearXZ     (atof(ARGUMENT));
    else if (OPTION("-syz")) dof.PutShearYZ     (atof(ARGUMENT));
    // Common pre-alignments
    else if (OPTION("-align-centers") || OPTION("-align-centres")) {
      align_centers = true;
    }
    else if (OPTION("-align-cog") || OPTION("-align-cogs")) {
      align_cogs = true;
    }
    else if (OPTION("-align-top")) {
      align_top = true;
    }
    // Input displacements (unstructured)
    else if (OPTION("-approximate") || OPTION("-approx")) {
      PointSetPair p;
      p.target_name = ARGUMENT;
      psets.push_back(p);
    }
    else if (OPTION("-displacements") || OPTION("-disp")) {
      PointSetPair p;
      p.target_name = ARGUMENT;
      if (HAS_ARGUMENT) {
        p.source_name = ARGUMENT;
        if (HAS_ARGUMENT) p.corrin_name = ARGUMENT;
      }
      psets.push_back(p);
    }
    // FFD lattice attributes
    else if (OPTION("-dx"))  sx = atof(ARGUMENT);
    else if (OPTION("-dy"))  sy = atof(ARGUMENT);
    else if (OPTION("-dz"))  sz = atof(ARGUMENT);
    else if (OPTION("-ds"))  sx = sy = sz = atof(ARGUMENT);
    else if (OPTION("-subdiv")) nsubdiv = atoi(ARGUMENT);
    // Target/source image
    else if (OPTION("-target")) target_name = ARGUMENT;
    else if (OPTION("-source")) source_name = ARGUMENT;
    else if (OPTION("-Tp")) target_padding = atof(ARGUMENT);
    else if (OPTION("-Sp")) source_padding = atof(ARGUMENT);
    // Common or unknown option
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (nsubdiv <  0) nsubdiv = 4;
  if (nsubdiv == 0) nsubdiv = 1;

  if (nonrigid && !target_name) {
    cerr << "Target image required when creating a non-rigid FFD!" << endl;
    exit(1);
  }

  // Read attributes of target image
  RealImage target;
  if (target_name) target.Read(target_name);
  const ImageAttributes &attr = target.Attributes();

  // ---------------------------------------------------------------------------
  // Usage 1) Common preset alignments
  if (align_top || align_centers || align_cogs) {

    if (!target_name || !source_name) {
      cerr << "Target and source image required for -align-* options" << endl;
      exit(1);
    }
    if (!psets.empty()) {
      cerr << "Warning: Ignoring input point sets when using -align-* options" << endl;
    }
    RealImage source(source_name);
    const ImageAttributes &sattr = source.Attributes();

    if (align_centers || align_cogs) {
      Matrix T(4, 4);
      T.Ident();
      Point tcenter, scenter;
      if (!align_cogs || IsNaN(target_padding) || target.CenterOfForeground(tcenter, target_padding) == 0) {
        tcenter = Point(attr._xorigin, attr._yorigin, attr._zorigin);
      }
      if (!align_cogs || IsNaN(source_padding) || source.CenterOfForeground(scenter, source_padding) == 0) {
        scenter = Point(sattr._xorigin, sattr._yorigin, sattr._zorigin);
      }
      dof.Transform(tcenter);
      T(0, 3) = scenter._x - tcenter._x;
      T(1, 3) = scenter._y - tcenter._y;
      T(2, 3) = scenter._z - tcenter._z;
      dof.PutMatrix(T * dof.GetMatrix());
    }

    if (align_top) {
      // Target center point of top slice
      Point tcenter, scenter;
      if (IsNaN(target_padding) || target.CenterOfForeground(tcenter, target_padding) == 0) {
        tcenter = Point(attr._xorigin, attr._yorigin, attr._zorigin);
      }
      target.WorldToImage(tcenter);
      tcenter._z = target.Z() - .5;
      target.ImageToWorld(tcenter);
      dof.Transform(tcenter);
      // Source center point of top slice
      if (IsNaN(source_padding) || source.CenterOfForeground(scenter, source_padding) == 0) {
        scenter = Point(sattr._xorigin, sattr._yorigin, sattr._zorigin);
      }
      target.WorldToImage(tcenter);
      tcenter._z = target.Z() - .5;
      target.ImageToWorld(tcenter);
      // Compose with translation to align top centers
      Matrix T(4, 4);
      T.Ident();
      T(0, 3) = scenter._x - tcenter._x;
      T(1, 3) = scenter._y - tcenter._y;
      T(2, 3) = scenter._z - tcenter._z;
      dof.PutMatrix(T * dof.GetMatrix());
    }

    if (verbose) dof.Print();
    dof.Write(dofout_name);

  }

  // ---------------------------------------------------------------------------
  // Usage 2) Create affine (FFD) transformation from parameter values read
  //          from command-line or input transformation file
  else if (psets.empty()) {

    // Initialize affine transformation from target and source image attributes
    if (source_name && dof.IsIdentity()) {
      GreyImage source(source_name);
      const ImageAttributes &sattr = source.Attributes();
      Point tcenter, scenter;
      if (IsNaN(target_padding) || target.CenterOfForeground(tcenter, target_padding) == 0) {
        tcenter = Point(attr._xorigin, attr._yorigin, attr._zorigin);
      }
      if (IsNaN(source_padding) || source.CenterOfForeground(scenter, source_padding) == 0) {
        scenter = Point(sattr._xorigin, sattr._yorigin, sattr._zorigin);
      }
      dof.PutTranslationX(scenter._x - tcenter._x);
      dof.PutTranslationY(scenter._y - tcenter._y);
      dof.PutTranslationZ(scenter._z - tcenter._z);
    }

    // Reset passive parameters
    AffineTransformation copy(dof);
    dof.CopyFrom(&copy);
    if (similarity) {
      dof.PutScale((dof.GetScaleX() + dof.GetScaleY() + dof.GetScaleZ()) / 3.0);
    }

    // Write affine transformation as global transformation of MFFD
    if (multilevel) {

      MultiLevelFreeFormTransformation mffd(dof);
      if (verbose) mffd.Print();
      mffd.Write(dofout_name);

    // Write affine transformation as free-form deformation
    } else if (nonrigid) {

      UniquePtr<FreeFormTransformation3D> ffd;
      if (diffeomorphic) ffd.reset(BSplineSVFFD(dof, attr, sx, sy, sz));
      else               ffd.reset(BSplineFFD  (dof, attr, sx, sy, sz));
      if (verbose) {
        dof.Print();
        ffd->Print();
      }
      ffd->Write(dofout_name);

    // Write affine transformation
    } else {

      if (verbose) dof.Print();
      dof.Write(dofout_name);

    }

  // ---------------------------------------------------------------------------
  // Usage 3) Create transformation which approximates given point displacements
  } else {

    // -------------------------------------------------------------------------
    // Combine input point sets into single unstructured displacement map

    PointSet target_points, source_points;
    HomogeneousTransformation *lin;
    FreeFormTransformation    *ffd;
    MultiLevelTransformation  *mffd;

    int n = 0;
    for (size_t i = 0; i < psets.size(); ++i) {
      if (psets[i].source_name == NULL && IsTransformation(psets[i].target_name)) {
        Transformation *dof = Transformation::New(psets[i].target_name);
        lin  = dynamic_cast<HomogeneousTransformation *>(dof);
        ffd  = dynamic_cast<FreeFormTransformation    *>(dof);
        mffd = dynamic_cast<MultiLevelTransformation  *>(dof);
        if (mffd && mffd->NumberOfLevels() > 0) ffd = mffd->GetLocalTransformation(-1);
        if (lin) {
          n += number_of_point_samples;
        } else if (ffd) {
          n += ffd->NumberOfActiveCPs();
        } else {
          cerr << "Unsupported transformation file: " << psets[i].target_name << endl;
          delete dof;
          exit(1);
        }
        delete dof;
      } else {
        target_points.Read(psets[i].target_name);
        n += target_points.Size();
      }
    }

    // Allocate memory for unstructured displacement map
    double *x  = Allocate<double>(n);
    double *y  = Allocate<double>(n);
    double *z  = Allocate<double>(n);
    double *dx = Allocate<double>(n);
    double *dy = Allocate<double>(n);
    double *dz = Allocate<double>(n);

    int j = 0;
    for (size_t i = 0; i < psets.size(); ++i) {

      // Add displacements of corresponding landmark points
      if (psets[i].source_name) {

        target_points.Read(psets[i].target_name);
        source_points.Read(psets[i].source_name);

        // Calculate unstructured displacement map
        if (psets[i].corrin_name) {
          Array<int> corr = ReadCorrespondences(psets[i].corrin_name,
                                                target_points.Size(),
                                                source_points.Size());
          for (int t = 0; t < target_points.Size(); ++t, ++j) {
            const int &s = corr[t];
            x [j] = target_points(t)._x;
            y [j] = target_points(t)._y;
            z [j] = target_points(t)._z;
            dx[j] = source_points(s)._x - x[j];
            dy[j] = source_points(s)._y - y[j];
            dz[j] = source_points(s)._z - z[j];
          }
        } else {
          if (source_points.Size() < target_points.Size()) {
            cerr << "Second point set has different number of points! An explicit correspondence map is required." << endl;
            Deallocate(x);
            Deallocate(y);
            Deallocate(z);
            Deallocate(dx);
            Deallocate(dy);
            Deallocate(dz);
            exit(1);
          }
          for (int t = 0; t < target_points.Size(); ++t, ++j) {
            x [j] = target_points(t)._x;
            y [j] = target_points(t)._y;
            z [j] = target_points(t)._z;
            dx[j] = source_points(t)._x - x[j];
            dy[j] = source_points(t)._y - y[j];
            dz[j] = source_points(t)._z - z[j];
          }
        }

      // Add point displacements from input transformation file
      } else if (IsTransformation(psets[i].target_name)) {

        Transformation *dof = Transformation::New(psets[i].target_name);
        lin  = dynamic_cast<HomogeneousTransformation *>(dof);
        ffd  = dynamic_cast<FreeFormTransformation    *>(dof);
        mffd = dynamic_cast<MultiLevelTransformation  *>(dof);
        if (mffd && mffd->NumberOfLevels() > 0) ffd = mffd->GetLocalTransformation(-1);
        if (lin) {

          PointSamples samples(number_of_point_samples);
          samples.SampleGaussian(.0, 100.0);

          int j1 = j;
          for (int i = 0; i < samples.Size(); ++i, ++j) {
            x[j] = samples(i)._x;
            y[j] = samples(i)._y;
            z[j] = samples(i)._z;
          }

          lin->Transform(samples);

          j = j1;
          for (int i = 0; i < samples.Size(); ++i, ++j) {
            dx[j] = samples(i)._x - x[j];
            dy[j] = samples(i)._y - y[j];
            dz[j] = samples(i)._z - z[j];
          }

        } else if (ffd) {

          GenericImage<double> disp(ffd->Attributes(), 3);
          ffd->Displacement(disp);

          for (int c = 0; c < disp.Z(); ++c)
          for (int b = 0; b < disp.Y(); ++b)
          for (int a = 0; a < disp.X(); ++a, ++j) {
            x[j] = a, y[j] = b, z[j] = c;
            disp.ImageToWorld(x[j], y[j], z[j]);
            dx[j] = disp(a, b, c, 0);
            dy[j] = disp(a, b, c, 1);
            dz[j] = disp(a, b, c, 2);
          }

        }
        delete dof;

      // Add point displacements read from VTK point set file
      } else {

        #if defined(HAVE_VTK) && MIRTK_IO_WITH_VTK
          vtkSmartPointer<vtkPointSet> pset = ReadPointSet(psets[i].target_name);
          bool pset_ok = true;
          if (pset->GetNumberOfPoints() == 0) {
            cerr << "Failed to read point set " << psets[i].target_name << " or it contains no points" << endl;
            pset_ok = false;
          }
          vtkDataArray *disp = pset->GetPointData()->GetVectors();
          if (pset_ok && disp == NULL) {
            cerr << "Point set " << psets[i].target_name << " has no displacement vectors!" << endl;
            pset_ok = false;
          }
          if (pset_ok && disp->GetNumberOfComponents() != 3) {
            cerr << "Point set " << psets[i].target_name << " displacement vectors must have dimension 3!" << endl;
            pset_ok = false;
          }
          if (!pset_ok) {
            Deallocate(x);
            Deallocate(y);
            Deallocate(z);
            Deallocate(dx);
            Deallocate(dy);
            Deallocate(dz);
            exit(1);
          }
          double p[3], v[3];
          for (vtkIdType ptId = 0; ptId < pset->GetNumberOfPoints(); ++ptId, ++j) {
            pset->GetPoint(ptId, p);
            disp->GetTuple(ptId, v);
            x [j] = p[0];
            y [j] = p[1];
            z [j] = p[2];
            dx[j] = v[0];
            dy[j] = v[1];
            dz[j] = v[2];
          }
        #else // MIRTK_IO_WITH_VTK
          cerr << EXECNAME << ": Link binary to PointSet module to be able to read displacements from VTK file!" << endl;
          exit(1);
        #endif // MIRTK_IO_WITH_VTK
      }
    }

    // -------------------------------------------------------------------------
    // 2a) Find non-rigid (M)FFD which approximates the input displacements
    if (nonrigid) {

      sx *= pow(2.0, nsubdiv-1);
      sy *= pow(2.0, nsubdiv-1);
      sz *= pow(2.0, nsubdiv-1);

      double *rx = Allocate<double>(n);
      double *ry = Allocate<double>(n);
      double *rz = Allocate<double>(n);

      if (multilevel) {

        memcpy(rx, dx, n * sizeof(double));
        memcpy(ry, dy, n * sizeof(double));
        memcpy(rz, dz, n * sizeof(double));
        dof.ApproximateAsNew(x, y, z, rx, ry, rz, n);

        MultiLevelFreeFormTransformation mffd(dof);
        FreeFormTransformation *ffd;
        if (diffeomorphic) ffd = new BSplineFreeFormTransformationSV(attr, sx, sy, sz);
        else               ffd = new BSplineFreeFormTransformation3D(attr, sx, sy, sz);
        mffd.PushLocalTransformation(ffd);

        double rms;
        for (int l = 0; l < nsubdiv; ++l) {
          if (l > 0) ffd->Subdivide();
          memcpy(rx, dx, n * sizeof(double));
          memcpy(ry, dy, n * sizeof(double));
          memcpy(rz, dz, n * sizeof(double));
          rms = ffd->Approximate(x, y, z, rx, ry, rz, n);
        }

        if (verbose) {
          cout << "RMS error of approximation is " << rms << "\n" << endl;
          mffd.Print();
        }
        mffd.Write(dofout_name);

      } else {

        UniquePtr<FreeFormTransformation3D> ffd;
        if (diffeomorphic) ffd.reset(new BSplineFreeFormTransformationSV(attr, sx, sy, sz));
        else               ffd.reset(new BSplineFreeFormTransformation3D(attr, sx, sy, sz));

        double rms;
        for (int l = 0; l < nsubdiv; ++l) {
          if (l > 0) ffd->Subdivide();
          memcpy(rx, dx, n * sizeof(double));
          memcpy(ry, dy, n * sizeof(double));
          memcpy(rz, dz, n * sizeof(double));
          rms = ffd->Approximate(x, y, z, rx, ry, rz, n);
        }

        if (verbose) {
          cout << "RMS error of approximation is " << rms << "\n" << endl;
          ffd->Print();
        }
        ffd->Write(dofout_name);

      }

      Deallocate(rx);
      Deallocate(ry);
      Deallocate(rz);

    // -------------------------------------------------------------------------
    // 2b) Find affine transformation which approximates the input displacements
    } else {

      // Minimize mean squared error of approximation
      double rms = dof.ApproximateAsNew(x, y, z, dx, dy, dz, n);
      if (verbose) cout << "RMS error of approximation is " << rms << "\n" << endl;

      if (similarity) {
        dof.PutScale((dof.GetScaleX() + dof.GetScaleY() + dof.GetScaleZ()) / 3.0);
      }

      // Write transformation
      if (multilevel) {
        MultiLevelFreeFormTransformation mffd(dof);
        mffd.Print();
        mffd.Write(dofout_name);
      } else {
        if (verbose) dof.Print();
        dof.Write(dofout_name);
      }

    }

    // -------------------------------------------------------------------------
    // Free displacement map
    Deallocate(x);
    Deallocate(y);
    Deallocate(z);
    Deallocate(dx);
    Deallocate(dy);
    Deallocate(dz);
  }

  return 0;
}
