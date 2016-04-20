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

// TODO: Rewrite test using GTest library. -as12312

#include "mirtk/GenericImage.h"
#include "mirtk/InterpolateImageFunction.h"
#include "mirtk/ExtrapolateImageFunction.h"

using namespace mirtk;

// ===========================================================================
// Macros
// ===========================================================================

// ---------------------------------------------------------------------------
#define TEST(id) \
  const char *_id    = id; \
  int         _nfail = 0; \
  cout << "Test case: " << _id << endl

// ---------------------------------------------------------------------------
#define RESULT _nfail

// ---------------------------------------------------------------------------
#define NULLPTR(actual, message, exit_if_not) \
  do { \
    void *_actual = (actual); \
    if (_actual != NULL) { \
      _nfail++; \
      cerr << message << endl; \
      cerr << "  Actual:   " << _actual << endl; \
      cerr << "  Expected: NULL" << endl; \
      if (exit_if_not) exit(_nfail); \
    } \
  }while(false)

// ---------------------------------------------------------------------------
#define NOT_NULLPTR(actual, message, exit_if_not) \
  do { \
    void *_actual = (actual); \
    if (_actual == NULL) { \
      _nfail++; \
      cerr << message << endl; \
      cerr << "  Actual:   " << _actual << endl; \
      cerr << "  Expected: !NULL" << endl; \
      if (exit_if_not) exit(_nfail); \
    } \
  }while(false)

// ---------------------------------------------------------------------------
#define EQUAL(actual, expected, message, exit_if_not) \
  do { \
    double _actual   = (actual); \
    double _expected = (expected); \
    if (fabs(_actual - _expected) >= 1e-10) { \
      _nfail++; \
      cerr << message << endl; \
      cerr << "  Actual:   " << setprecision(15) << _actual   << endl; \
      cerr << "  Expected: " << setprecision(15) << _expected << endl; \
      if (exit_if_not) exit(_nfail); \
    } \
  }while(false)

// ---------------------------------------------------------------------------
#define NOT_EQUAL(actual, expected, message, exit_if_not) \
  do { \
    double _actual   = (actual); \
    double _expected = (expected); \
    if (fabs(_actual - _expected) < 1e-10) { \
      _nfail++; \
      cerr << message << endl; \
      cerr << "  Actual:   " << setprecision(15) << _actual   << endl; \
      cerr << "  Expected: " << setprecision(15) << _expected << endl; \
      if (exit_if_not) exit(_nfail); \
    } \
  }while(false)

// ---------------------------------------------------------------------------
#define STREQUAL(actual, expected, message, exit_if_not) \
  do { \
    string _actual   = (actual); \
    string _expected = (expected); \
    if (_actual != _expected) { \
      _nfail++; \
      cerr << message << endl; \
      cerr << "  Actual:   " << _actual << endl; \
      cerr << "  Expected: " << _expected << endl; \
      if (exit_if_not) exit(_nfail); \
    } \
  }while(false)

// ---------------------------------------------------------------------------
#define EXPECT_NULL(actual, message)                NULLPTR(actual, message, false)
#define EXPECT_NOT_NULL(actual, message)            NOT_NULLPTR(actual, message, false)
#define EXPECT_EQUAL(actual, expected, message)     EQUAL(actual, expected, message, false)
#define EXPECT_NOT_EQUAL(actual, expected, message) NOT_EQUAL(actual, expected, message, false)
#define EXPECT_STREQUAL (actual, expected, message) STREQUAL(actual, expected, message, false)

// ---------------------------------------------------------------------------
#define ASSERT_NULL(actual, message)                NULLPTR(actual, message, true)
#define ASSERT_NOT_NULL(actual, message)            NOT_NULLPTR(actual, message, true)
#define ASSERT_EQUAL(actual, expected, message)     EQUAL(actual, expected, message, true)
#define ASSERT_NOT_EQUAL(actual, expected, message) NOT_EQUAL(actual, expected, message, true)
#define ASSERT_STREQUAL(actual, expected, message)  STREQUAL(actual, expected, message, true)

// ===========================================================================
// Test images
// ===========================================================================

// ---------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType> *create_test_image(int x = 16, int y = 16, int z = 1, int t = 1, int n = 1)
{
  // Instantiate image
  GenericImage<VoxelType> *image = new GenericImage<double>(x, y, z, t, n);
  // Fill test image
  for (int l = 0; l < image->T(); l++) {
    for (int k = 0; k < image->Z(); k++) {
      for (int j = 0; j < image->Y(); j++) {
        for (int i = 0; i < image->X(); i++) {
          image->Put(i, j, k, l, static_cast<VoxelType>(image->VoxelToIndex(i, j, k, l)));
        }
      }
    }
  }
  return image;
}

// ===========================================================================
// Tests
// ===========================================================================

// ---------------------------------------------------------------------------
int test_Initialize()
{
  TEST("test_Initialize");

  GenericImage<double>     *image = create_test_image<double>(16, 16);
  InterpolateImageFunction *inter = InterpolateImageFunction::New(Interpolation_NN, Extrapolation_Const, image);
  ExtrapolateImageFunction *extra = inter->Extrapolator();

  inter->Initialize();
  extra->Initialize(); // not required, but shall also not hurt (ie., recurse infinitely)
  inter->Initialize();

  return RESULT;
}

// ---------------------------------------------------------------------------
int test_IsInsideOutside()
{
  TEST("test_IsInsideOutside");

  GenericImage<double>     *image = NULL;
  InterpolateImageFunction *inter = NULL;
  //ExtrapolateImageFunction *extra = NULL;

  // --------------------------------------------------------------------------
  // 2D
  image = create_test_image<double>(16, 16);
  inter = InterpolateImageFunction::New(Interpolation_NN, Extrapolation_Const, image);
  //extra = inter->Extrapolator();

  inter->Initialize();

  EXPECT_EQUAL(inter->IsInside (3.0, 2.0), true,  "Check if voxel is inside image domain");
  EXPECT_EQUAL(inter->IsInside (3.4, 2.1), true,  "Check if sub-voxel is inside image domain");

  EXPECT_EQUAL(inter->IsInside (-3.0, 0.0), false, "Check if voxel is outside image domain");
  EXPECT_EQUAL(inter->IsInside (-3.4, 2.1), false, "Check if sub-voxel is outside image domain");

  delete image;
  delete inter;

  // --------------------------------------------------------------------------
  // 3D
  image = create_test_image<double>(16, 16, 16);
  inter = InterpolateImageFunction::New(Interpolation_NN, Extrapolation_Const, image);
  //extra = inter->Extrapolator();

  inter->Initialize();

  EXPECT_EQUAL(inter->IsInside (3.0, 2.0, 5.0), true,  "Check if voxel is inside image domain");
  EXPECT_EQUAL(inter->IsInside (3.4, 2.1, 5.5), true,  "Check if sub-voxel is inside image domain");

  EXPECT_EQUAL(inter->IsInside (-3.0, 0.0, 1.0), false, "Check if voxel is outside image domain");
  EXPECT_EQUAL(inter->IsInside (-3.4, 2.1, 1.0), false, "Check if sub-voxel is outside image domain");

  delete image;
  delete inter;

  return RESULT;
}

// ---------------------------------------------------------------------------
int test_Interpolation_NN_Extrapolation_Default()
{
  TEST("test_Interpolation_NN_Extrapolation_Default");

  GenericImage<double>     *image = create_test_image<double>(16, 16, 16);
  InterpolateImageFunction *func  = InterpolateImageFunction::New(Interpolation_NN, image);

  // Initialize image function
  EXPECT_NULL(func->Extrapolator(), "Extrapolator not instantiated yet");

  func->Initialize();

  ASSERT_NOT_NULL(func->Extrapolator(), "Extrapolator initialized");

  // Check instantiated image function
  ASSERT_STREQUAL(func->NameOfClass(), "GenericNearestNeighborInterpolateImageFunction", "Type of interpolator");
  ASSERT_STREQUAL(func->Extrapolator()->NameOfClass(), "GenericConstExtrapolateImageFunction", "Type of extrapolator");
  EXPECT_EQUAL(func->Extrapolator()->DefaultValue(), 0, "Constant outside value");

  // Interpolate inside of image domain
  EXPECT_EQUAL(func->Evaluate(3.0, 2.0, 5.0), image->VoxelToIndex(3, 2, 5), "Interpolate image value at voxel position");  
  EXPECT_EQUAL(func->Evaluate(3.4, 2.1, 5.5), image->VoxelToIndex(3, 2, 6), "Interpolate image value at sub-voxel position");

  // Extrapolate outside of image domain
  EXPECT_EQUAL(func->Evaluate(-3.0, 0.0, 1.0), 0, "Extrapolate image value at voxel position");  
  EXPECT_EQUAL(func->Evaluate(-3.4, 2.1, 1.0), 0, "Extrapolate image value at sub-voxel position");

  // Clean up
  delete image;
  delete func;

  return RESULT;
}

// ---------------------------------------------------------------------------
int test_Interpolation_NN_Extrapolation_Const()
{
  TEST("test_Interpolation_NN_Extrapolation_Const");

  GenericImage<double>     *image = NULL;
  InterpolateImageFunction *func  = NULL;

  // -------------------------------------------------------------------------
  // 2D
  image = create_test_image<double>(16, 16);
  func  = InterpolateImageFunction::New(Interpolation_NN, Extrapolation_Const, image);

  // Initialize image function
  func->Initialize();

  // Check instantiated image function
  ASSERT_STREQUAL(func->NameOfClass(), "GenericNearestNeighborInterpolateImageFunction", "Type of interpolator");
  ASSERT_STREQUAL(func->Extrapolator()->NameOfClass(), "GenericConstExtrapolateImageFunction", "Type of extrapolator");
  EXPECT_EQUAL(func->Extrapolator()->DefaultValue(), 0, "Constant outside value");

  // Interpolate inside of image domain
  EXPECT_EQUAL(func->Evaluate(3.0, 2.0), image->VoxelToIndex(3, 2), "Interpolate image value at voxel position");  
  EXPECT_EQUAL(func->Evaluate(3.4, 2.1), image->VoxelToIndex(3, 2), "Interpolate image value at sub-voxel position");

  // Extrapolate outside of image domain
  EXPECT_EQUAL(func->Evaluate(-3.0, 0.0), 0, "Extrapolate image value at voxel position");  
  EXPECT_EQUAL(func->Evaluate(-3.4, 2.1), 0, "Extrapolate image value at sub-voxel position");

  // Clean up
  delete image;
  delete func;
  
  // -------------------------------------------------------------------------
  // 3D
  image = create_test_image<double>(16, 16, 16);
  func  = InterpolateImageFunction::New(Interpolation_NN, Extrapolation_Const, image);

  // Initialize image function
  func->Initialize();

  // Check instantiated image function
  ASSERT_STREQUAL(func->NameOfClass(), "GenericNearestNeighborInterpolateImageFunction", "Type of interpolator");
  ASSERT_STREQUAL(func->Extrapolator()->NameOfClass(), "GenericConstExtrapolateImageFunction", "Type of extrapolator");
  EXPECT_EQUAL(func->Extrapolator()->DefaultValue(), 0, "Constant outside value");

  // Interpolate inside of image domain
  EXPECT_EQUAL(func->Evaluate(3.0, 2.0, 5.0), image->VoxelToIndex(3, 2, 5), "Interpolate image value at voxel position");  
  EXPECT_EQUAL(func->Evaluate(3.4, 2.1, 5.5), image->VoxelToIndex(3, 2, 6), "Interpolate image value at sub-voxel position");

  // Extrapolate outside of image domain
  EXPECT_EQUAL(func->Evaluate(-3.0, 0.0, 1.0), 0, "Extrapolate image value at voxel position");  
  EXPECT_EQUAL(func->Evaluate(-3.4, 2.1, 1.0), 0, "Extrapolate image value at sub-voxel position");

  // Clean up
  delete image;
  delete func;

  return RESULT;
}

// ---------------------------------------------------------------------------
int test_Interpolation_Linear_Extrapolation_Const()
{
  TEST("test_Interpolation_Linear_Extrapolation_Const");

  GenericImage<double>     *image = NULL;
  InterpolateImageFunction *func  = NULL;
  
  // -------------------------------------------------------------------------
  // 2D
  image = create_test_image<double>(16, 16);
  func  = InterpolateImageFunction::New(Interpolation_Linear, Extrapolation_Const, image);

  // Initialize image function
  func->Initialize();
  func->Extrapolator()->DefaultValue(-1);

  // Check instantiated image function
  ASSERT_STREQUAL(func->NameOfClass(), "GenericLinearInterpolateImageFunction2D", "Type of interpolator");
  ASSERT_STREQUAL(func->Extrapolator()->NameOfClass(), "GenericConstExtrapolateImageFunction", "Type of extrapolator");

  // Interpolate inside image domain
  EXPECT_EQUAL(func->Evaluate(3.0, 2.0), image->VoxelToIndex(3, 2), "Interpolate image value at voxel position");  
  EXPECT_EQUAL(func->Evaluate(3.4, 2.0), 0.6 * image->VoxelToIndex(3, 2) + 0.4 * image->VoxelToIndex(4, 2), "Interpolate image value at sub-voxel position along x");
  EXPECT_EQUAL(func->Evaluate(3.0, 2.7), 0.3 * image->VoxelToIndex(3, 2) + 0.7 * image->VoxelToIndex(3, 3), "Interpolate image value at sub-voxel position along y");

  // Extrapolate outside image domain
  ASSERT_EQUAL(func->Extrapolator()->DefaultValue(), -1, "Default value of extrapolator");
  EXPECT_EQUAL(func->Evaluate(-0.1, -0.1), -0.19, "Extrapolate image value outside of image domain");
  EXPECT_EQUAL(func->Evaluate(-1.1,  0.0), -1, "Extrapolate image value outside of image domain in x");
  EXPECT_EQUAL(func->Evaluate(16.1,  0.0), -1, "Extrapolate image value outside of image domain in x");
  EXPECT_EQUAL(func->Evaluate( 0.0, -1.1), -1, "Extrapolate image value outside of image domain in y");
  EXPECT_EQUAL(func->Evaluate( 0.0, 16.1), -1, "Extrapolate image value outside of image domain in y");
  EXPECT_EQUAL(func->Evaluate(-1.1, -1.1), -1, "Extrapolate image value outside of image domain in x and y");
  EXPECT_EQUAL(func->Evaluate(16.1, 16.1), -1, "Extrapolate image value outside of image domain in x and y");

  // Clean up
  delete image;
  delete func;

  // -------------------------------------------------------------------------
  // 3D
  image = create_test_image<double>(16, 16, 16);
  func  = InterpolateImageFunction::New(Interpolation_Linear, Extrapolation_Const, image);

  // Initialize image function
  func->Initialize();
  func->Extrapolator()->DefaultValue(-1);

  // Check instantiated image function
  ASSERT_STREQUAL(func->NameOfClass(), "GenericLinearInterpolateImageFunction3D", "Type of interpolator");
  ASSERT_STREQUAL(func->Extrapolator()->NameOfClass(), "GenericConstExtrapolateImageFunction", "Type of extrapolator");

  // Interpolate inside image domain
  EXPECT_EQUAL(func->Evaluate(3.0, 2.0, 5.0), image->VoxelToIndex(3, 2, 5), "Interpolate image value at voxel position");  
  EXPECT_EQUAL(func->Evaluate(3.4, 2.0, 5.0), 0.6 * image->VoxelToIndex(3, 2, 5) + 0.4 * image->VoxelToIndex(4, 2, 5), "Interpolate image value at sub-voxel position along x");
  EXPECT_EQUAL(func->Evaluate(3.0, 2.7, 5.0), 0.3 * image->VoxelToIndex(3, 2, 5) + 0.7 * image->VoxelToIndex(3, 3, 5), "Interpolate image value at sub-voxel position along y");
  EXPECT_EQUAL(func->Evaluate(3.0, 2.0, 5.1), 0.9 * image->VoxelToIndex(3, 2, 5) + 0.1 * image->VoxelToIndex(3, 2, 6), "Interpolate image value at sub-voxel position along z");

  // Extrapolate outside image domain
  ASSERT_EQUAL(func->Extrapolator()->DefaultValue(), -1, "Default value of extrapolator");
  EXPECT_EQUAL(func->Evaluate(-1.1,  0.0,  0.0), -1, "Extrapolate image value slightly outside of image domain in x");
  EXPECT_EQUAL(func->Evaluate(16.1,  0.0,  0.0), -1, "Extrapolate image value slightly outside of image domain in x");
  EXPECT_EQUAL(func->Evaluate( 0.0, -1.1,  0.0), -1, "Extrapolate image value slightly outside of image domain in y");
  EXPECT_EQUAL(func->Evaluate( 0.0, 16.1,  0.0), -1, "Extrapolate image value slightly outside of image domain in y");
  EXPECT_EQUAL(func->Evaluate( 0.0,  0.0, -1.1), -1, "Extrapolate image value slightly outside of image domain in z");
  EXPECT_EQUAL(func->Evaluate( 0.0,  0.0, 16.1), -1, "Extrapolate image value slightly outside of image domain in z");
  EXPECT_EQUAL(func->Evaluate(-1.1, -1.1,  0.0), -1, "Extrapolate image value slightly outside of image domain in x and y");
  EXPECT_EQUAL(func->Evaluate(16.1, 16.1,  0.0), -1, "Extrapolate image value slightly outside of image domain in x and y");
  EXPECT_EQUAL(func->Evaluate(-1.1,  0.0, -1.1), -1, "Extrapolate image value slightly outside of image domain in x and z");
  EXPECT_EQUAL(func->Evaluate(16.1,  0.0, 16.1), -1, "Extrapolate image value slightly outside of image domain in x and z");
  EXPECT_EQUAL(func->Evaluate( 0.0, -1.1, -1.1), -1, "Extrapolate image value slightly outside of image domain in y and z");
  EXPECT_EQUAL(func->Evaluate( 0.0, 16.1, 16.1), -1, "Extrapolate image value slightly outside of image domain in y and z");

  // Clean up
  delete image;
  delete func;

  // -------------------------------------------------------------------------
  // 3D vector field
  image = create_test_image<double>(16, 16, 16, 1, 3);
  func  = InterpolateImageFunction::New(Interpolation_Linear, Extrapolation_Const, image);

  // Initialize image function
  func->Initialize();
  func->Extrapolator()->DefaultValue(-1);

  // Check instantiated image function
  ASSERT_STREQUAL(func->NameOfClass(), "GenericLinearInterpolateImageFunction3D", "Type of interpolator");
  ASSERT_STREQUAL(func->Extrapolator()->NameOfClass(), "GenericConstExtrapolateImageFunction", "Type of extrapolator");

  // Interpolate inside image domain
  EXPECT_EQUAL(func->Evaluate(3.0, 2.0, 5.0, 0.0), image->VoxelToIndex(3, 2, 5, 0), "Interpolate image value at voxel position");  
  EXPECT_EQUAL(func->Evaluate(3.4, 2.0, 5.0, 0.0), 0.6 * image->VoxelToIndex(3, 2, 5, 0) + 0.4 * image->VoxelToIndex(4, 2, 5, 0), "Interpolate image value at sub-voxel position along x");
  EXPECT_EQUAL(func->Evaluate(3.0, 2.7, 5.0, 0.0), 0.3 * image->VoxelToIndex(3, 2, 5, 0) + 0.7 * image->VoxelToIndex(3, 3, 5, 0), "Interpolate image value at sub-voxel position along y");
  EXPECT_EQUAL(func->Evaluate(3.0, 2.0, 5.1, 0.0), 0.9 * image->VoxelToIndex(3, 2, 5, 0) + 0.1 * image->VoxelToIndex(3, 2, 6, 0), "Interpolate image value at sub-voxel position along z");

  EXPECT_EQUAL(func->Evaluate(3.0, 2.0, 5.0, 1.0), image->VoxelToIndex(3, 2, 5, 1), "Interpolate image value at voxel position");  
  EXPECT_EQUAL(func->Evaluate(3.4, 2.0, 5.0, 1.0), 0.6 * image->VoxelToIndex(3, 2, 5, 1) + 0.4 * image->VoxelToIndex(4, 2, 5, 1), "Interpolate image value at sub-voxel position along x");
  EXPECT_EQUAL(func->Evaluate(3.0, 2.7, 5.0, 1.0), 0.3 * image->VoxelToIndex(3, 2, 5, 1) + 0.7 * image->VoxelToIndex(3, 3, 5, 1), "Interpolate image value at sub-voxel position along y");
  EXPECT_EQUAL(func->Evaluate(3.0, 2.0, 5.1, 1.0), 0.9 * image->VoxelToIndex(3, 2, 5, 1) + 0.1 * image->VoxelToIndex(3, 2, 6, 1), "Interpolate image value at sub-voxel position along z");

  EXPECT_EQUAL(func->Evaluate(3.0, 2.0, 5.0, 2.0), image->VoxelToIndex(3, 2, 5, 2), "Interpolate image value at voxel position");  
  EXPECT_EQUAL(func->Evaluate(3.4, 2.0, 5.0, 2.0), 0.6 * image->VoxelToIndex(3, 2, 5, 2) + 0.4 * image->VoxelToIndex(4, 2, 5, 2), "Interpolate image value at sub-voxel position along x");
  EXPECT_EQUAL(func->Evaluate(3.0, 2.7, 5.0, 2.0), 0.3 * image->VoxelToIndex(3, 2, 5, 2) + 0.7 * image->VoxelToIndex(3, 3, 5, 2), "Interpolate image value at sub-voxel position along y");
  EXPECT_EQUAL(func->Evaluate(3.0, 2.0, 5.1, 2.0), 0.9 * image->VoxelToIndex(3, 2, 5, 2) + 0.1 * image->VoxelToIndex(3, 2, 6, 2), "Interpolate image value at sub-voxel position along z");

  // Clean up
  delete image;
  delete func;

  return RESULT;
}

// ---------------------------------------------------------------------------
int test_Interpolation_Linear_Extrapolation_NN()
{
  TEST("test_Interpolation_Linear_Extrapolation_NN");

  GenericImage<double>     *image = NULL;
  InterpolateImageFunction *func  = NULL;
  
  // -------------------------------------------------------------------------
  // 2D
  image = create_test_image<double>(16, 16);
  func  = InterpolateImageFunction::New(Interpolation_Linear, Extrapolation_NN, image);

  // Initialize image function
  func->Initialize();

  // Check instantiated image function
  ASSERT_STREQUAL(func->NameOfClass(), "GenericLinearInterpolateImageFunction2D", "Type of interpolator");
  ASSERT_STREQUAL(func->Extrapolator()->NameOfClass(), "GenericNearestNeighborExtrapolateImageFunction", "Type of extrapolator");

  // Interpolate inside image domain
  EXPECT_EQUAL(func->Evaluate(3.0, 2.0), image->VoxelToIndex(3, 2), "Interpolate image value at voxel position");  
  EXPECT_EQUAL(func->Evaluate(3.4, 2.0), 0.6 * image->VoxelToIndex(3, 2) + 0.4 * image->VoxelToIndex(4, 2), "Interpolate image value at sub-voxel position along x");
  EXPECT_EQUAL(func->Evaluate(3.0, 2.7), 0.3 * image->VoxelToIndex(3, 2) + 0.7 * image->VoxelToIndex(3, 3), "Interpolate image value at sub-voxel position along y");

  // Extrapolate outside image domain
  EXPECT_EQUAL(func->Evaluate(-0.1,  0.0), image->VoxelToIndex( 0,  0), "Extrapolate image value slightly outside of image domain in x");
  EXPECT_EQUAL(func->Evaluate(15.1,  0.0), image->VoxelToIndex(15,  0), "Extrapolate image value slightly outside of image domain in x");
  EXPECT_EQUAL(func->Evaluate( 0.0, -0.1), image->VoxelToIndex( 0,  0), "Extrapolate image value slightly outside of image domain in y");
  EXPECT_EQUAL(func->Evaluate( 0.0, 15.1), image->VoxelToIndex( 0, 15), "Extrapolate image value slightly outside of image domain in y");
  EXPECT_EQUAL(func->Evaluate(-0.1, -0.1), image->VoxelToIndex( 0,  0), "Extrapolate image value slightly outside of image domain in x and y");
  EXPECT_EQUAL(func->Evaluate(15.1, 15.1), image->VoxelToIndex(15, 15), "Extrapolate image value slightly outside of image domain in x and y");

  // Clean up
  delete image;
  delete func;

  // -------------------------------------------------------------------------
  // 3D
  image = create_test_image<double>(16, 16, 16);
  func  = InterpolateImageFunction::New(Interpolation_Linear, Extrapolation_NN, image);

  // Initialize image function
  func->Initialize();

  // Check instantiated image function
  ASSERT_STREQUAL(func->NameOfClass(), "GenericLinearInterpolateImageFunction3D", "Type of interpolator");
  ASSERT_STREQUAL(func->Extrapolator()->NameOfClass(), "GenericNearestNeighborExtrapolateImageFunction", "Type of extrapolator");

  // Interpolate inside image domain
  EXPECT_EQUAL(func->Evaluate(3.0, 2.0, 5.0), image->VoxelToIndex(3, 2, 5), "Interpolate image value at voxel position");  
  EXPECT_EQUAL(func->Evaluate(3.4, 2.0, 5.0), 0.6 * image->VoxelToIndex(3, 2, 5) + 0.4 * image->VoxelToIndex(4, 2, 5), "Interpolate image value at sub-voxel position along x");
  EXPECT_EQUAL(func->Evaluate(3.0, 2.7, 5.0), 0.3 * image->VoxelToIndex(3, 2, 5) + 0.7 * image->VoxelToIndex(3, 3, 5), "Interpolate image value at sub-voxel position along y");
  EXPECT_EQUAL(func->Evaluate(3.0, 2.0, 5.1), 0.9 * image->VoxelToIndex(3, 2, 5) + 0.1 * image->VoxelToIndex(3, 2, 6), "Interpolate image value at sub-voxel position along z");

  // Extrapolate outside image domain
  EXPECT_EQUAL(func->Evaluate(-0.1,  0.0,  0.0), image->VoxelToIndex( 0,  0,  0), "Extrapolate image value slightly outside of image domain in x");
  EXPECT_EQUAL(func->Evaluate(15.1,  0.0,  0.0), image->VoxelToIndex(15,  0,  0), "Extrapolate image value slightly outside of image domain in x");
  EXPECT_EQUAL(func->Evaluate( 0.0, -0.1,  0.0), image->VoxelToIndex( 0,  0,  0), "Extrapolate image value slightly outside of image domain in y");
  EXPECT_EQUAL(func->Evaluate( 0.0, 15.1,  0.0), image->VoxelToIndex( 0, 15,  0), "Extrapolate image value slightly outside of image domain in y");
  EXPECT_EQUAL(func->Evaluate( 0.0,  0.0, -0.1), image->VoxelToIndex( 0,  0,  0), "Extrapolate image value slightly outside of image domain in z");
  EXPECT_EQUAL(func->Evaluate( 0.0,  0.0, 15.1), image->VoxelToIndex( 0,  0, 15), "Extrapolate image value slightly outside of image domain in z");
  EXPECT_EQUAL(func->Evaluate(-0.1, -0.1,  0.0), image->VoxelToIndex( 0,  0,  0), "Extrapolate image value slightly outside of image domain in x and y");
  EXPECT_EQUAL(func->Evaluate(15.1, 15.1,  0.0), image->VoxelToIndex(15, 15,  0), "Extrapolate image value slightly outside of image domain in x and y");
  EXPECT_EQUAL(func->Evaluate(-0.1,  0.0, -0.1), image->VoxelToIndex( 0,  0,  0), "Extrapolate image value slightly outside of image domain in x and z");
  EXPECT_EQUAL(func->Evaluate(15.1,  0.0, 15.1), image->VoxelToIndex(15,  0, 15), "Extrapolate image value slightly outside of image domain in x and z");
  EXPECT_EQUAL(func->Evaluate( 0.0, -0.1, -0.1), image->VoxelToIndex( 0,  0,  0), "Extrapolate image value slightly outside of image domain in y and z");
  EXPECT_EQUAL(func->Evaluate( 0.0, 15.1, 15.1), image->VoxelToIndex( 0, 15, 15), "Extrapolate image value slightly outside of image domain in y and z");

  // Clean up
  delete image;
  delete func;

  return RESULT;
}

// ===========================================================================
// Main
// ===========================================================================

// ---------------------------------------------------------------------------
int InterpolateExtrapolateImageFunctionTest(int, char *[])
{
  int retval = 0;

  retval += test_Initialize();
  retval += test_IsInsideOutside();
  retval += test_Interpolation_NN_Extrapolation_Default();
  retval += test_Interpolation_NN_Extrapolation_Const();
  retval += test_Interpolation_Linear_Extrapolation_Const();
  retval += test_Interpolation_Linear_Extrapolation_NN();

  return retval;
}
