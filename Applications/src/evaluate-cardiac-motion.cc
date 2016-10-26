/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Wenjia Bai
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
#include "mirtk/Transformations.h"
#include "mirtk/PointSetIO.h"

#include "mirtk/Vtk.h"
#include "mirtk/VtkMath.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataNormals.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkCenterOfMass.h"


using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
// Print help screen
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <input_mesh> <output_mesh> -image <file> -dof <file>" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Evaluate cardiac motion. This command computes the displacement or strain" << endl;
  cout << "  at each vertex on a myocardial surface mesh." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input_mesh            Myocardial mesh." << endl;
  cout << "  output_mesh           The output mesh which stores the vertex-wise motion data." << endl;
  cout << endl;
  cout << "Required options:" << endl;
  cout << "  -image <file>         Cardiac image whose z-axis defines the longitudinal direction." << endl;
  cout << "  -dof, -dofin <file>   Transformation with regard to the reference frame (end-diastolic frame)." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -invert-zaxis         Use this option if the longitudinal direction is the inverted z-axis of the image. (default: off)" << endl;
  cout << "  -displacement         Compute the displacement. (default: off)" << endl;
  cout << "  -strain               Compute the strain.       (default: off)" << endl;
  cout << "  -save-local-coord     Save the local cardiac coordinate system in the output file. (default: off)" << endl;
  cout << "  -ascii                Write legacy VTK files encoded in ASCII. (default: off)" << endl;
  cout << "  -binary               Write legacy VTK files in binary form. (default: on)" << endl;
  cout << "  -compress             Compress XML VTK files. (default: on)" << endl;
  cout << "  -nocompress           Do not compress XML VTK files. (default: off)" << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Determine local coordinate system for strain calculations
///
/// - Caroline Petitjean et al. Assessment of myocardial function: a review of
///   quantification methods and results using tagged MRI.
///   Journal of Cardiovascular Magnetic Resonance, 2005.
/// - An Elen et al. Three-Dimensional Cardiac Strain Estimation Using Spatio-Temporal
///   Elastic Registration of Ultrasound Images: A Feasibility Study.
///   IEEE Transactions on Medical Imaging, 2008.
class CalculateLocalDirections
{
  Point         _Center;
  vtkPoints    *_Points;
  vtkDataArray *_Normals;
  vtkDataArray *_Radial;
  vtkDataArray *_Longit;
  vtkDataArray *_Circum;

public:

  /// Operator called be parallel_for, not to be called directly!
  void operator ()(const blocked_range<vtkIdType> &ptIds) const
  {
    double cs, dp, p[3], v_pc[3], radial[3], longit[3], circum[3];;
    for (vtkIdType ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      // Get point coordinates
      _Points->GetPoint(ptId, p);

      // The longitudinal direction is defined to be pointing from apex to base.
      // Normally, we assume that z-axis of the input image (short-axis image stack) is the longitudinal axis.
      _Longit->GetTuple(ptId, longit);
   
      // The radial direction is defined to be pointing outside the epicardium.
      // The normal vectors for the endocardium computed by the vtkPolyDataNormals filter point outside, which is ok for us.
      // However, the normal vectors for the epicardium may point inside, which need to be flipped.
      // Check whether normal_vec is pointing outside or not. If it points inside, invert it.
      // v_pc = c - p
      _Normals->GetTuple(ptId, radial);
      v_pc[0] = _Center._x - p[0];
      v_pc[1] = _Center._y - p[1];
      v_pc[2] = _Center._z - p[2];
      
      // Remove the longitudinal component of v_pc
      dp = vtkMath::Dot(v_pc, longit);
      v_pc[0] -= dp * longit[0];
      v_pc[1] -= dp * longit[1];
      v_pc[2] -= dp * longit[2];

      if (vtkMath::Dot(radial, v_pc) > 0.) {
        vtkMath::MultiplyScalar(radial, -1.);
      }

      // At the apex, the cardiac local coordinate system is illy defined, since radial can
      // be in the same direction as longit. Also, the circumferential direction is also undefined at the apex.
      // In this case, we set both radial and circumferential direction to be zero.
      cs = abs(vtkMath::Dot(radial, longit)) / (vtkMath::Norm(radial) * vtkMath::Norm(longit));
      if (fequal(cs, 1.)) {
        // Extreme case at the apex
        radial[0] = circum[0] = 0.;
        radial[1] = circum[1] = 0.;
        radial[2] = circum[2] = 0.;
      } else {
        // Normal case
        // Make sure the radial direction is perpendicular to the longitudinal direction
        // Subtract its projection on the longitudinal direction and normalise it to a unit vector.
        dp = vtkMath::Dot(radial, longit);
        radial[0] -= dp * longit[0];
        radial[1] -= dp * longit[1];
        radial[2] -= dp * longit[2];
        vtkMath::Normalize(radial);

        // Determine the circumferential direction which is the cross product of the radial and the longitudinal directions
        vtkMath::Cross(radial, longit, circum);
      }
      
      // Store the local coordinate system
      _Radial->SetTuple(ptId, radial);
      _Circum->SetTuple(ptId, circum);
    }
  }

  /// Calculate local coordinate axes for each mesh node
  static void Run(vtkPolyData           *mesh,
                  const ImageAttributes &attr,
                  vtkDataArray          *radial,
                  vtkDataArray          *longit,
                  vtkDataArray          *circum)
  {
    CalculateLocalDirections body;

    // Compute the centre of left ventricle
    double c[3];
    vtkNew<vtkCenterOfMass> center_of_mass;
    SetVTKInput(center_of_mass, mesh);
    center_of_mass->SetUseScalarsAsWeights(false);
    center_of_mass->Update();
    center_of_mass->GetCenter(c);
    body._Center = Point(c);

    // Compute surface normals
    vtkSmartPointer<vtkDataArray> normals;
    vtkNew<vtkPolyDataNormals> normals_filter;
    normals_filter->SetInputData(mesh);
    normals_filter->SetFeatureAngle(30);
    normals_filter->SetSplitting(0);         
    normals_filter->SetConsistency(1);
    normals_filter->SetAutoOrientNormals(1);
    normals_filter->SetFlipNormals(0);
    normals_filter->SetNonManifoldTraversal(1);
    normals_filter->SetComputePointNormals(1);
    normals_filter->SetComputeCellNormals(0);
    normals_filter->Update();
    normals = normals_filter->GetOutput()->GetPointData()->GetNormals();
    body._Normals = normals;

    // The longitudinal direction is defined to be pointing from apex to base.
    // Normally, we assume that z-axis of the input image (short-axis image stack) is the longitudinal axis.
    longit->FillComponent(0, attr._zaxis[0]);
    longit->FillComponent(1, attr._zaxis[1]);
    longit->FillComponent(2, attr._zaxis[2]);

    // Compute other directions
    body._Points = mesh->GetPoints();
    body._Radial = radial;
    body._Longit = longit;
    body._Circum = circum;
    parallel_for(blocked_range<vtkIdType>(0, mesh->GetNumberOfPoints()), body);
  }
};

// -----------------------------------------------------------------------------
/// Functor for computing transformed points
struct TransformPoints
{
  vtkPoints            *_Points;
  const Transformation *_Transform;
  vtkPoints            *_TransformedPoints;

  void operator ()(const blocked_range<vtkIdType> ptIds) const
  {
    double p[3];
    for (vtkIdType ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      _Points->GetPoint(ptId, p);
      _Transform->Transform(p[0], p[1], p[2]);
      _TransformedPoints->SetPoint(ptId, p);
    }
  }
};

// -----------------------------------------------------------------------------
/// Functor for computation of vertex displacements
struct CalculateDirectionalDisplacements
{
  vtkPoints    *_Points;
  vtkPoints    *_TransformedPoints;
  vtkDataArray *_Radial;
  vtkDataArray *_Longit;
  vtkDataArray *_Circum;
  vtkDataArray *_Displacement;

  void operator ()(const blocked_range<vtkIdType> ptIds) const
  {
    double p[3], q[3], d[3], radial[3], longit[3], circum[3];
    for (vtkIdType ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {

      // Get point coordinates before and after transformation
      _Points->GetPoint(ptId, p);
      _TransformedPoints->GetPoint(ptId, q);

      // Get directions of local coordinate system
      _Radial->GetTuple(ptId, radial);
      _Longit->GetTuple(ptId, longit);
      _Circum->GetTuple(ptId, circum);

      // Calculate local displacement vectors
      vtkMath::Subtract(q, p, d);
      _Displacement->SetComponent(ptId, 0, vtkMath::Dot(d, radial));
      _Displacement->SetComponent(ptId, 1, vtkMath::Dot(d, longit));
      _Displacement->SetComponent(ptId, 2, vtkMath::Dot(d, circum));
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate strain along the local coordinate axes
///
/// Transform the Jacobian matrix to the local coordinate system, compute tensor
/// matrix, extract normal strains: \f$E = 0.5 * (F^T \dot F - I)\f$.
///
/// - Caroline Petitjean et al. Assessment of myocardial function: a review of quantification
///   methods and results using tagged MRI. Journal of Cardiovascular Magnetic Resonance, 2005.
///   Equation (3)
/// - Lewis K. Waldman et al. Transmural myocardial deformation in the canine left ventricle.
///   Circulation Research, Vol. 57, No. 1, July 1985. Equation (8)
struct CalculateDirectionalStrain
{
  vtkPoints            *_Points;
  const Transformation *_Transform;
  vtkDataArray         *_Radial;
  vtkDataArray         *_Longit;
  vtkDataArray         *_Circum;
  vtkDataArray         *_Strain;

  void operator ()(const blocked_range<vtkIdType> ptIds) const
  {
    double p[3], radial[3], longit[3], circum[3];
    Matrix J(3, 3), J_l(3, 3), P(3, 3), E_l(3, 3);
    for (vtkIdType ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {

      // Get point coordinates at which to evaluate Jacobian of transformation
      _Points->GetPoint(ptId, p);

      // Get directions of local coordinate system
      _Radial->GetTuple(ptId, radial);
      _Longit->GetTuple(ptId, longit);
      _Circum->GetTuple(ptId, circum);
      P(0, 0) = radial[0]; P(0, 1) = longit[0]; P(0, 2) = circum[0];
      P(1, 0) = radial[1]; P(1, 1) = longit[1]; P(1, 2) = circum[1];
      P(2, 0) = radial[2]; P(2, 1) = longit[2]; P(2, 2) = circum[2];

      // Calculate Jacobian of transformation in local coordinate system
      _Transform->Jacobian(J, p[0], p[1], p[2]);
      J_l  = P.Transposed();
      J_l *= J;
      J_l *= P;

      // Calculate strain tensor
      E_l  = J_l.Transposed();
      E_l *= J_l;

      E_l(0, 0) -= 1.;
      E_l(1, 1) -= 1.;
      E_l(2, 2) -= 1.;

      E_l *= 0.5;

      // Normal strains in three directions, which represent the change of squared difference.
      // Shear strains are not reported.
      _Strain->SetComponent(ptId, 0, E_l(0, 0));
      _Strain->SetComponent(ptId, 1, E_l(1, 1));
      _Strain->SetComponent(ptId, 2, E_l(2, 2));
    }
  }
};

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  EXPECTS_POSARGS(2);

  const char *input_mesh_name  = POSARG(1);
  const char *output_mesh_name = POSARG(2);
  FileOption  output_mesh_fopt = FO_Default;
  const char *image_name       = nullptr;
  const char *dof_name         = nullptr;

  bool invert_zaxis     = false;
  bool compute_disp     = false;
  bool compute_strain   = false;
  bool save_local_coord = false;

  for (ALL_OPTIONS) {
    if      (OPTION("-image")) image_name = ARGUMENT;
    else if (OPTION("-dof") || OPTION("-dofin")) dof_name = ARGUMENT;
    else if (OPTION("-invert-zaxis")) invert_zaxis = true;
    else if (OPTION("-displacement")) compute_disp = true;
    else if (OPTION("-strain")) compute_strain = true;
    else if (OPTION("-save-local-coord")) save_local_coord = true;
    else HANDLE_POINTSETIO_OPTION(output_mesh_fopt);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (image_name == nullptr) {
    FatalError("Input -image file name required!");
  }
  if (dof_name == nullptr) {
    FatalError("Input -dof file name required!");
  }

  // Read input image
  ImageAttributes attr;
  {
    if (verbose) cout << "Reading image attributes from " << image_name << "...";
    InitializeIOLibrary();
    GreyImage image(image_name);
    attr = image.Attributes();
    if (invert_zaxis) {
      // Normally, we assume that z-axis of the input image (short-axis image stack)
      // is the longitudinal axis, but sometimes this assumption may be wrong
      // and we need to invert the z-axis
      attr._zaxis[0] = -attr._zaxis[0];
      attr._zaxis[1] = -attr._zaxis[1];
      attr._zaxis[2] = -attr._zaxis[2];
    }
    if (verbose) cout << " done" << endl;
  }

  // Read input mesh
  if (verbose) cout << "Reading mesh from " << input_mesh_name << "...";
  FileOption input_mesh_fopt;
  vtkSmartPointer<vtkPolyData> input_mesh;
  input_mesh = ReadPolyData(input_mesh_name, input_mesh_fopt);
  if (verbose) cout << " done" << endl;

  const vtkIdType npoints = input_mesh->GetNumberOfPoints();
  blocked_range<vtkIdType> ptIds(0, npoints);

  // Read input transformation
  if (verbose) cout << "Reading transformation from " << dof_name << "...";
  UniquePtr<Transformation> dof(Transformation::New(dof_name));
  if (verbose) cout << " done" << endl;

  // Compute deformed output mesh
  vtkSmartPointer<vtkPolyData> output_mesh;
  output_mesh.TakeReference(input_mesh->NewInstance());
  output_mesh->ShallowCopy(input_mesh);
  vtkSmartPointer<vtkPoints> output_points;
  output_points = vtkSmartPointer<vtkPoints>::New();
  output_points->SetNumberOfPoints(npoints);
  output_mesh->SetPoints(output_points);
  if (output_mesh_fopt == FO_Default) output_mesh_fopt = input_mesh_fopt;

  TransformPoints transform;
  transform._Points            = input_mesh->GetPoints();
  transform._Transform         = dof.get();
  transform._TransformedPoints = output_mesh->GetPoints();
  parallel_for(ptIds, transform);

  // Determine the local cardiac coordinate system for each point
  vtkSmartPointer<vtkDataArray> radial, longit, circum;

  radial = vtkSmartPointer<vtkFloatArray>::New();
  radial->SetName("Dir_Radial");
  radial->SetNumberOfComponents(3);
  radial->SetNumberOfTuples(npoints);

  longit = vtkSmartPointer<vtkFloatArray>::New();
  longit->SetName("Dir_Longit");
  longit->SetNumberOfComponents(3);
  longit->SetNumberOfTuples(npoints);

  circum = vtkSmartPointer<vtkFloatArray>::New();
  circum->SetName("Dir_Circum");
  circum->SetNumberOfComponents(3);
  circum->SetNumberOfTuples(npoints);

  // Use the input mesh (the reference mesh) to determine the local coordinate system
  CalculateLocalDirections::Run(input_mesh, attr, radial, longit, circum);

  if (save_local_coord) {
    output_mesh->GetPointData()->AddArray(radial);
    output_mesh->GetPointData()->AddArray(longit);
    output_mesh->GetPointData()->AddArray(circum);
  }

  // Calculate local displacements
  if (compute_disp) {
    vtkSmartPointer<vtkDataArray> disp;
    disp = vtkSmartPointer<vtkFloatArray>::New();
    disp->SetName("Displacement");
    disp->SetNumberOfComponents(3);
    disp->SetNumberOfTuples(npoints);
    output_mesh->GetPointData()->AddArray(disp);

    CalculateDirectionalDisplacements calc_disp;
    calc_disp._Points            = input_mesh->GetPoints();
    calc_disp._TransformedPoints = output_mesh->GetPoints();
    calc_disp._Radial            = radial;
    calc_disp._Longit            = longit;
    calc_disp._Circum            = circum;
    calc_disp._Displacement      = disp;
    parallel_for(ptIds, calc_disp);
  }

  // Calculate local strain
  if (compute_strain) {
    vtkSmartPointer<vtkDataArray> strain;
    strain = vtkSmartPointer<vtkFloatArray>::New();
    strain->SetName("Strain");
    strain->SetNumberOfComponents(3);
    strain->SetNumberOfTuples(npoints);
    output_mesh->GetPointData()->AddArray(strain);

    CalculateDirectionalStrain calc_strain;
    calc_strain._Points     = input_mesh->GetPoints();
    calc_strain._Transform  = dof.get();
    calc_strain._Radial     = radial;
    calc_strain._Longit     = longit;
    calc_strain._Circum     = circum;
    calc_strain._Strain     = strain;
    parallel_for(ptIds, calc_strain);
  }

  // Write the output mesh
  if (!WritePolyData(output_mesh_name, output_mesh, output_mesh_fopt)) {
    FatalError("Failed to write output mesh to file " << output_mesh_name);
  }

  return 0;
}
