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

#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/IOConfig.h"
#include "mirtk/GenericImage.h"
#include "mirtk/Transformations.h"

#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataNormals.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkCenterOfMass.h"
#include "vtkPolyDataWriter.h"


using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
// Print help screen
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <image> <mesh> <T> <mesh_out>" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Evaluate cardiac motion. This command computes the displacement or strain" << endl;
  cout << "  at each vertex on a myocardial surface mesh." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  image:    input image whose z-axis provides the longitudinal direction." << endl;
  cout << "  mesh:     input myocardial mesh." << endl;
  cout << "  T:        transformation with regard to the reference frame (end-diastolic frame)." << endl;
  cout << "  mesh_out: output mesh which stores the vertex-wise motion data." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -invert-zaxis      Use this option if the longitudinal direction is the inverted z-axis of the image. (Default: off)" << endl;
  cout << "  -displacement      Compute the displacement. (Default: off)" << endl;
  cout << "  -strain            Compute the strain.       (Default: off)" << endl;
  cout << "  -save-local-coord  Save the local cardiac coordinate system in the output file. (Default: off)" << endl;
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


// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  REQUIRES_POSARGS(4);
  const char *image_name    = POSARG(1);
  const char *mesh_name     = POSARG(2);
  const char *T_name        = POSARG(3);
  const char *mesh_out_name = POSARG(4);
  bool  invert_zaxis        = false;
  bool  compute_disp        = false;
  bool  compute_strain      = false;
  bool  save_local_coord    = false;
  
  for (ALL_OPTIONS) {
    if (OPTION("-invert-zaxis")) {
      invert_zaxis = true;
    }
    else if (OPTION("-displacement")) {
      compute_disp = true;
    }
    else if (OPTION("-strain")) {
      compute_strain = true;
    }
    else if (OPTION("-save-local-coord")) {
      save_local_coord = true;
    }
    else {
      HANDLE_STANDARD_OR_UNKNOWN_OPTION();
    }
  }

  // Read input image
  InitializeIOLibrary();
  GreyImage image(image_name);
  if (verbose) { cout << "Reading image from " << image_name << endl; }

  ImageAttributes attr = image.Attributes();

  // Read input mesh
  vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName(mesh_name);
  reader->Update();
  if (verbose) { cout << "Reading mesh from " << mesh_name << endl; }

  vtkSmartPointer<vtkPolyData> mesh = reader->GetOutput();
  vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
  int n_points = points->GetNumberOfPoints();

  // Read input transformation
  Transformation *T = dynamic_cast<Transformation *>(Transformation::New(T_name));
  if (verbose) { cout << "Reading transformation from " << T_name << endl; }

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
  vtkSmartPointer<vtkFloatArray> disp = vtkSmartPointer<vtkFloatArray>::New();
  disp->SetName("Displacement");
  disp->SetNumberOfComponents(3);
  disp->SetNumberOfTuples(n_points);

  vtkSmartPointer<vtkFloatArray> strain = vtkSmartPointer<vtkFloatArray>::New();
  strain->SetName("Strain");
  strain->SetNumberOfComponents(3);
  strain->SetNumberOfTuples(n_points);

  // For each point on the mesh
  for (int i = 0; i < n_points; i++) {
    // Point coordinate
    double p[3];
    points->GetPoint(i, p);

    // Get the local coordinate system at this point
    double radial[3], longit[3], circum[3];
    dir_radial->GetTuple(i, radial);
    dir_longit->GetTuple(i, longit);
    dir_circum->GetTuple(i, circum);

    // Transform this point
    double q[3] = {p[0], p[1], p[2]};
    T->Transform(q[0], q[1], q[2]);

    // Update the point coordinate
    points->SetPoint(i, q);

    // Compute displacement
    double d[3];
    for (int dim = 0; dim < 3; dim++) {
      d[dim] = q[dim] - p[dim];
    }

    // Displacement in local coordinate system
    double d_l[3] = {vtkMath::Dot(d, radial),
		     vtkMath::Dot(d, longit),
		     vtkMath::Dot(d, circum)};
    disp->SetTuple(i, d_l);

    // Compute strain
    // Method: transform the Jacobian matrix to the local coordinate system, compute tensor matrix, extract normal strains.
    Matrix J;
    T->Jacobian(J, p[0], p[1], p[2]);

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
    strain->SetComponent(i, 0, E_l(0,0));
    strain->SetComponent(i, 1, E_l(1,1));
    strain->SetComponent(i, 2, E_l(2,2));
  }

  // Save the displacement and strain data
  if (compute_disp)        { mesh->GetPointData()->AddArray(disp); }
  if (compute_strain)      { mesh->GetPointData()->AddArray(strain); }
  if (save_local_coord) {
    mesh->GetPointData()->AddArray(dir_radial);
    mesh->GetPointData()->AddArray(dir_longit);
    mesh->GetPointData()->AddArray(dir_circum);
  }
  
  // Write the output mesh
  if (verbose) { cout << "Writing mesh to " << mesh_out_name << endl; }
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInputData(mesh);
  writer->SetFileName(mesh_out_name);
  writer->Update();

  return 0;
}
