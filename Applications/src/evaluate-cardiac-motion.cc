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
#include "mirtk/VtkMath.h"

#include "vtkDataReader.h"
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
  cout << "  input_mesh           Myocardial mesh." << endl;
  cout << "  output_mesh          The output mesh which stores the vertex-wise motion data." << endl;
  cout << "  -image file          Cardiac image whose z-axis defines the longitudinal direction." << endl;
  cout << "  -dof, -dofin file    Transformation with regard to the reference frame (end-diastolic frame)." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -invert-zaxis        Use this option if the longitudinal direction is the inverted z-axis of the image. (default: off)" << endl;
  cout << "  -displacement        Compute the displacement. (default: off)" << endl;
  cout << "  -strain              Compute the strain.       (default: off)" << endl;
  cout << "  -save-local-coord    Save the local cardiac coordinate system in the output file. (default: off)" << endl; 
  cout << "  -ascii               Write legacy VTK files encoded in ASCII. (default: off)" << endl;
  cout << "  -binary              Write legacy VTK files in binary form. (default: on)" << endl;
  cout << "  -compress            Compress XML VTK files. (default: on)" << endl;
  cout << "  -nocompress          Do not compress XML VTK files. (default: off)" << endl;
 PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
// Determine the cardiac local coordinate system
//
// References:
// [1] Caroline Petitjean et al. Assessment of myocardial function: a review of quantification methods and results using tagged MRI. Journal of Cardiovascular Magnetic Resonance, 2005. 
// [2] An Elen et al. Three-Dimensional Cardiac Strain Estimation Using Spatio-Temporal Elastic Registration of Ultrasound Images: A Feasibility Study. IEEE Transactions on Medical Imaging, 2008.
void DetermineLocalCoordinateSystem(vtkSmartPointer<vtkPolyData> mesh, ImageAttributes attr, bool invert_zaxis,
                                    vtkSmartPointer<vtkDataArray> dir_radial,
                                    vtkSmartPointer<vtkDataArray> dir_longit,
                                    vtkSmartPointer<vtkDataArray> dir_circum)
{
  // Prepare the array to store the local coordinate system
  vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
  int n_points = points->GetNumberOfPoints();

  // Compute the centre of left ventricle
  vtkSmartPointer<vtkCenterOfMass> center_of_mass = vtkSmartPointer<vtkCenterOfMass>::New();
  double c[3];
  center_of_mass->SetInputData(mesh);
  center_of_mass->SetUseScalarsAsWeights(0);
  center_of_mass->Update();
  center_of_mass->GetCenter(c);

  // Compute surface normals
  vtkSmartPointer<vtkPolyDataNormals> normals_filter = vtkSmartPointer<vtkPolyDataNormals>::New();
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
  vtkSmartPointer<vtkDataArray> normals = normals_filter->GetOutput()->GetPointData()->GetNormals();

  // For each point on the mesh
  for (int i = 0; i < n_points; i++) {
    // Point coordinate
    double p[3];
    points->GetPoint(i, p);

    // Determine the local coordinate system at this point
    double radial[3], longit[3], circum[3];

    // The longitudinal direction is defined to be pointing from apex to base.
    // Normally, we assume that z-axis of the input image (short-axis image stack) is the longitudinal axis.
    // NOTE, sometimes this assumption may be wrong and we need to invert the z-axis.
    if (invert_zaxis) {
      for (int dim = 0; dim < 3; dim++) {
        longit[dim] = - attr._zaxis[dim];
      }
    }
    else {
      for (int dim = 0; dim < 3; dim++) {
        longit[dim] = attr._zaxis[dim];
      }
    }
 
    // The radial direction is defined to be pointing outside the epicardium.
    // The normal vectors for the endocardium computed by the vtkPolyDataNormals filter point outside, which is ok for us.
    // However, the normal vectors for the epicardium may point inside, which need to be flipped.
    double normal_vec[3];
    normals->GetTuple(i, normal_vec);

    // Check whether normal_vec is pointing outside or not. If it points inside, invert it.
    // v_pc = c - p
    double v_pc[3];
    vtkMath::Subtract(c, p, v_pc);
    if (vtkMath::Dot(normal_vec, v_pc) > 0) {
      vtkMath::MultiplyScalar(normal_vec, -1);
    }

    // The radial direction is the outward normal vector
    for (int dim = 0; dim < 3; dim++) {
      radial[dim] = normal_vec[dim];
    }
      
    // At the apex, the cardiac local coordinate system is illy defined, since radial can be in the same direction as longit. Also, the circumferential direction is also undefined at the apex.
    // In this case, we set both radial and circumferential direction to be zero.
    double cos = fabs(vtkMath::Dot(radial, longit)) / (vtkMath::Norm(radial) * vtkMath::Norm(longit));
    if (cos == 1) {
      // Extreme case at the apex
      for (int dim = 0; dim < 3; dim++) {
        radial[dim] = 0;
        circum[dim] = 0;
      }
    }
    else {
      // Normal case
      // Make sure the radial direction is perpendicular to the longitudinal direction
      // Subtract its projection on the longitudinal direction and normalise it to a unit vector.
      double dot_prod = vtkMath::Dot(radial, longit);
      for (int dim = 0; dim < 3; dim++) {
        radial[dim] -= dot_prod * longit[dim];
      }
      vtkMath::Normalize(radial);

      // Determine the circumferential direction which is the cross product of the radial and the longitudinal directions
      vtkMath::Cross(radial, longit, circum);
    }
    
    // Store the local coordinate system
    dir_radial->SetTuple(i, radial);
    dir_longit->SetTuple(i, longit);
    dir_circum->SetTuple(i, circum);
  }
}

// -----------------------------------------------------------------------------
// Strain evaluation class
struct CalculateDisplacementAndStrain
{
  vtkPoints            *_Points;
  const Transformation *_Transform;
  vtkDataArray         *_Dir_Radial;
  vtkDataArray         *_Dir_Longit;
  vtkDataArray         *_Dir_Circum;
  vtkDataArray         *_Displacement;
  vtkDataArray         *_Strain;
  vtkPoints            *_TransformedPoints;

  void operator ()(const blocked_range<vtkIdType> ptIds) const
  {
    // For each point
    for (vtkIdType i = ptIds.begin(); i != ptIds.end(); ++i) {
      // Get point coordinate
      double p[3];
      _Points->GetPoint(i, p);

      // Get directions of local coordinate system
      double radial[3], longit[3], circum[3];
      _Dir_Radial->GetTuple(i, radial);
      _Dir_Longit->GetTuple(i, longit);
      _Dir_Circum->GetTuple(i, circum);

      // Transform point
      double q[3] = {p[0], p[1], p[2]};
      _Transform->Transform(q[0], q[1], q[2]);

      if (_TransformedPoints) {
        _TransformedPoints->SetPoint(i, q);
      }

      // Calculate local displacement vectors
      if (_Displacement) {
        double d[3];
        vtkMath::Subtract(q, p, d);
        _Displacement->SetComponent(i, 0, vtkMath::Dot(d, radial));
        _Displacement->SetComponent(i, 1, vtkMath::Dot(d, longit));
        _Displacement->SetComponent(i, 2, vtkMath::Dot(d, circum));
      }

      // Calculate strain
      // Method: transform the Jacobian matrix to the local coordinate system, compute tensor matrix, extract normal strains.
      if (_Strain) {
        Matrix J(3,3);
        _Transform->Jacobian(J, p[0], p[1], p[2]);

        Matrix P(3,3);
        P(0,0) = radial[0]; P(0,1) = longit[0]; P(0,2) = circum[0]; 
        P(1,0) = radial[1]; P(1,1) = longit[1]; P(1,2) = circum[1]; 
        P(2,0) = radial[2]; P(2,1) = longit[2]; P(2,2) = circum[2]; 

        // Jacobian in local coordinate system
        Matrix J_l(3,3);

        J_l = P;
        J_l.Transpose();
        J_l *= J;
        J_l *= P;

        // Strain tensor in local coordinate system
        // Math:
        //   E = 0.5 * (F^T \dot F - I)
        //
        // References:
        // [1] Caroline Petitjean et al. Assessment of myocardial function: a review of quantification methods and results using tagged MRI. Journal of Cardiovascular Magnetic Resonance, 2005. 
        //   Equation (3)
        // or an even older paper
        // [2] Lewis K. Waldman et al. Transmural myocardial deformation in the canine left ventricle. Circulation Research, Vol. 57, No. 1, July 1985.
        //   Equation (8)
        //
        Matrix E_l(3,3);
        E_l = J_l;
        E_l.Transpose();
        E_l *= J_l;
      
        Matrix I(3,3);
        I.Ident();

        E_l -= I;
        E_l *= 0.5;
      
        // Normal strains in three directions, which represent the change of squared difference.
        // Shear strains are not reported.
        _Strain->SetComponent(i, 0, E_l(0,0));
        _Strain->SetComponent(i, 1, E_l(1,1));
        _Strain->SetComponent(i, 2, E_l(2,2));
      }
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
  const char *image_name       = nullptr;
  const char *dof_name         = nullptr;
  bool  invert_zaxis     = false;
  bool  compute_disp     = false;
  bool  compute_strain   = false;
  bool  save_local_coord = false;
  int   ascii            = -1;
  bool  compress         = true;

  for (ALL_OPTIONS) {
    if      (OPTION("-image"))            { image_name = ARGUMENT; }
    else if (OPTION("-dof") || OPTION("-dofin")) { dof_name = ARGUMENT; }
    else if (OPTION("-invert-zaxis"))     { invert_zaxis = true; }
    else if (OPTION("-displacement"))     { compute_disp = true; }
    else if (OPTION("-strain"))           { compute_strain = true; }
    else if (OPTION("-save-local-coord")) { save_local_coord = true; }
    else if (OPTION("-ascii") || OPTION("-nobinary")) { ascii = 1; }
    else if (OPTION("-noascii") || OPTION("-binary")) { ascii = 0; }
    else if (OPTION("-compress"))         { compress = true; }
    else if (OPTION("-nocompress"))       { compress = false; } 
    else { HANDLE_STANDARD_OR_UNKNOWN_OPTION(); }
  }

  // Read input image
  if (image_name == nullptr) { FatalError("Input -image file name required!"); }
  ImageAttributes attr;
  {
    if (verbose) { cout << "Reading image attributes from " << image_name << "..."; }
    InitializeIOLibrary();
    GreyImage image(image_name);
    attr = image.Attributes();
    if (verbose) { cout << " done" << endl; }
  }

  // Read input mesh
  if (verbose) { cout << "Reading mesh from " << input_mesh_name << "..."; }
  int mesh_ftype;
  vtkSmartPointer<vtkPolyData> mesh = ReadPolyData(input_mesh_name, &mesh_ftype);
  if (verbose) { cout << " done" << endl; }
  if (ascii == -1) { ascii = (mesh_ftype == VTK_ASCII ? 1 : 0); }

  vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
  int n_points = points->GetNumberOfPoints();

  // Read input transformation
  if (dof_name == nullptr) { FatalError("Input -dof file name required!"); }
  if (verbose) { cout << "Reading transformation from " << dof_name << "..."; }
  UniquePtr<Transformation> T(Transformation::New(dof_name));
  if (verbose) { cout << " done" << endl; }

  // Determine the cardiac local coordinate system
  //
  // References:
  // [1] Caroline Petitjean et al. Assessment of myocardial function: a review of quantification methods and results using tagged MRI. Journal of Cardiovascular Magnetic Resonance, 2005. 
  // [2] An Elen et al. Three-Dimensional Cardiac Strain Estimation Using Spatio-Temporal Elastic Registration of Ultrasound Images: A Feasibility Study. IEEE Transactions on Medical Imaging, 2008.
  //
  vtkSmartPointer<vtkDataArray> dir_radial = vtkSmartPointer<vtkFloatArray>::New();
  dir_radial->SetName("Dir_Radial");
  dir_radial->SetNumberOfComponents(3);
  dir_radial->SetNumberOfTuples(n_points);

  vtkSmartPointer<vtkDataArray> dir_longit = vtkSmartPointer<vtkFloatArray>::New();
  dir_longit->SetName("Dir_Longit");
  dir_longit->SetNumberOfComponents(3);
  dir_longit->SetNumberOfTuples(n_points);

  vtkSmartPointer<vtkDataArray> dir_circum = vtkSmartPointer<vtkFloatArray>::New();
  dir_circum->SetName("Dir_Circum");
  dir_circum->SetNumberOfComponents(3);
  dir_circum->SetNumberOfTuples(n_points);

  DetermineLocalCoordinateSystem(mesh, attr, invert_zaxis, dir_radial, dir_longit, dir_circum);

  // Prepare the array for storing the displacements and strains
  vtkSmartPointer<vtkDataArray> disp;
  if (compute_disp) {
    disp = vtkSmartPointer<vtkFloatArray>::New();
    disp->SetName("Displacement");
    disp->SetNumberOfComponents(3);
    disp->SetNumberOfTuples(n_points);
  }

  vtkSmartPointer<vtkDataArray> strain;
  if (compute_strain) {
    strain = vtkSmartPointer<vtkFloatArray>::New();
    strain->SetName("Strain");
    strain->SetNumberOfComponents(3);
    strain->SetNumberOfTuples(n_points);
  }

  // Prepare the evaluation object
  CalculateDisplacementAndStrain eval;
  eval._Points            = points;
  eval._Transform         = T.get();
  eval._Dir_Radial        = dir_radial;
  eval._Dir_Longit        = dir_longit;
  eval._Dir_Circum        = dir_circum;
  eval._Displacement      = disp;
  eval._Strain            = strain;
  eval._TransformedPoints = points;

  // Evaluate the displacement and strain
  parallel_for(blocked_range<vtkIdType>(0, n_points), eval);

  // Save the displacement and strain data
  if (compute_disp)     { mesh->GetPointData()->AddArray(disp); }
  if (compute_strain)   { mesh->GetPointData()->AddArray(strain); }
  if (save_local_coord) {
    mesh->GetPointData()->AddArray(dir_radial);
    mesh->GetPointData()->AddArray(dir_longit);
    mesh->GetPointData()->AddArray(dir_circum);
  }
  
  // Write the output mesh
  if (!WritePolyData(output_mesh_name, mesh, compress, static_cast<bool>(ascii))) {
    FatalError("Failed to write output mesh to file " << output_mesh_name);
  }

  return 0;
}
