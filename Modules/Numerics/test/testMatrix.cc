/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#include "NumericsTest.h"

#include "mirtk/Matrix.h"
using namespace mirtk;

// =============================================================================
// Test matrices
// =============================================================================

// -----------------------------------------------------------------------------
Matrix AffineTestMatrix()
{
  Matrix m(4, 4);
  m(0, 0) =  0.8259, m(0, 1) =  0.1532, m(0, 2) = -0.1097, m(0, 3) =  -1.0848;
  m(1, 0) = -0.0641, m(1, 1) =  0.8929, m(1, 2) =  0.4143, m(1, 3) = -29.0270;
  m(2, 0) =  0.1467, m(2, 1) = -0.3754, m(2, 2) =  0.8963, m(2, 3) =  -1.5679;
  m(3, 0) =  0.0000, m(3, 1) =  0.0000, m(3, 2) =  0.0000, m(3, 3) =   1.0000;
  return m;
}

// =============================================================================
// Maps
// =============================================================================

// -----------------------------------------------------------------------------
TEST(Matrix, Sqrt)
{
  Matrix a = AffineTestMatrix();
  Matrix b = a.Sqrt();
  Matrix c = b * b;
  EXPECT_MATRIX_EQ(a, c);
}

// -----------------------------------------------------------------------------
TEST(Matrix, LogExp)
{
  Matrix a = AffineTestMatrix();
  Matrix b = a.Log();
  Matrix c = b.Exp();
  EXPECT_MATRIX_EQ(a, c);
}

// -----------------------------------------------------------------------------
TEST(Matrix, ExpLog)
{
  Matrix a = AffineTestMatrix();
  Matrix b = a.Exp();
  Matrix c = b.Log();
  EXPECT_MATRIX_EQ(a, c);
}

// =============================================================================
// Transformation matrices
// =============================================================================

// -----------------------------------------------------------------------------
TEST(Matrix, IdentityToRigidParameters)
{
  Matrix m;
  double tx, ty, tz, rx, ry, rz;
  RigidParametersToMatrix(.0, .0, .0, .0, .0, .0, m);
  MatrixToRigidParameters(m, tx, ty, tz, rx, ry, rz);
  EXPECT_DOUBLE_EQ( 0.0, tx);
  EXPECT_DOUBLE_EQ( 0.0, ty);
  EXPECT_DOUBLE_EQ( 0.0, tz);
  EXPECT_DOUBLE_EQ( 0.0, rx);
  EXPECT_DOUBLE_EQ( 0.0, ry);
  EXPECT_DOUBLE_EQ( 0.0, rz);
}

// -----------------------------------------------------------------------------
TEST(Matrix, TranslationToRigidParameters)
{
  Matrix m;
  double tx, ty, tz, rx, ry, rz;
  RigidParametersToMatrix(-1.0, 0.5, 4.7, .0, .0, .0, m);
  MatrixToRigidParameters(m, tx, ty, tz, rx, ry, rz);
  EXPECT_DOUBLE_EQ(-1.0, tx);
  EXPECT_DOUBLE_EQ( 0.5, ty);
  EXPECT_DOUBLE_EQ( 4.7, tz);
  EXPECT_DOUBLE_EQ( 0.0, rx);
  EXPECT_DOUBLE_EQ( 0.0, ry);
  EXPECT_DOUBLE_EQ( 0.0, rz);
}

// -----------------------------------------------------------------------------
TEST(Matrix, RotationToRigidParameters)
{
  Matrix m;
  double tx, ty, tz, rx, ry, rz;
  RigidParametersToMatrix(.0, .0, .0, .7, -1.1, 3.0, m);
  MatrixToRigidParameters(m, tx, ty, tz, rx, ry, rz);
  EXPECT_DOUBLE_EQ( 0.0, tx);
  EXPECT_DOUBLE_EQ( 0.0, ty);
  EXPECT_DOUBLE_EQ( 0.0, tz);
  EXPECT_DOUBLE_EQ( 0.7, rx);
  EXPECT_DOUBLE_EQ(-1.1, ry);
  EXPECT_DOUBLE_EQ( 3.0, rz);
}

// -----------------------------------------------------------------------------
TEST(Matrix, TranslationAndRotationToRigidParameters)
{
  Matrix m;
  double tx, ty, tz, rx, ry, rz;
  RigidParametersToMatrix(-1.0, 0.5, 4.7, .7, -1.1, 3.0, m);
  MatrixToRigidParameters(m, tx, ty, tz, rx, ry, rz);
  EXPECT_DOUBLE_EQ(-1.0, tx);
  EXPECT_DOUBLE_EQ( 0.5, ty);
  EXPECT_DOUBLE_EQ( 4.7, tz);
  EXPECT_DOUBLE_EQ( 0.7, rx);
  EXPECT_DOUBLE_EQ(-1.1, ry);
  EXPECT_DOUBLE_EQ( 3.0, rz);
}

// -----------------------------------------------------------------------------
TEST(Matrix, IdentityToAffineParameters)
{
  Matrix m;
  double tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz;
  AffineParametersToMatrix(.0, .0, .0, .0, .0, .0, 1.0, 1.0, 1.0, .0, .0, .0, m);
  MatrixToAffineParameters(m, tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz);
  EXPECT_DOUBLE_EQ( 0.0, tx);
  EXPECT_DOUBLE_EQ( 0.0, ty);
  EXPECT_DOUBLE_EQ( 0.0, tz);
  EXPECT_DOUBLE_EQ( 0.0, rx);
  EXPECT_DOUBLE_EQ( 0.0, ry);
  EXPECT_DOUBLE_EQ( 0.0, rz);
  EXPECT_DOUBLE_EQ( 1.0, sx);
  EXPECT_DOUBLE_EQ( 1.0, sy);
  EXPECT_DOUBLE_EQ( 1.0, sz);
  EXPECT_DOUBLE_EQ( 0.0, sxy);
  EXPECT_DOUBLE_EQ( 0.0, sxz);
  EXPECT_DOUBLE_EQ( 0.0, syz);
}

// -----------------------------------------------------------------------------
TEST(Matrix, TranslationToAffineParameters)
{
  Matrix m;
  double tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz;
  AffineParametersToMatrix(-1.0, 0.5, 4.7, .0, .0, .0, 1.0, 1.0, 1.0, .0, .0, .0, m);
  MatrixToAffineParameters(m, tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz);
  EXPECT_DOUBLE_EQ(-1.0, tx);
  EXPECT_DOUBLE_EQ( 0.5, ty);
  EXPECT_DOUBLE_EQ( 4.7, tz);
  EXPECT_DOUBLE_EQ( 0.0, rx);
  EXPECT_DOUBLE_EQ( 0.0, ry);
  EXPECT_DOUBLE_EQ( 0.0, rz);
  EXPECT_DOUBLE_EQ( 1.0, sx);
  EXPECT_DOUBLE_EQ( 1.0, sy);
  EXPECT_DOUBLE_EQ( 1.0, sz);
  EXPECT_DOUBLE_EQ( 0.0, sxy);
  EXPECT_DOUBLE_EQ( 0.0, sxz);
  EXPECT_DOUBLE_EQ( 0.0, syz);
}

// -----------------------------------------------------------------------------
TEST(Matrix, RotationToAffineParameters)
{
  Matrix m;
  double tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz;
  AffineParametersToMatrix(.0, .0, .0, .7, -1.1, 3.0, 1.0, 1.0, 1.0, .0, .0, .0, m);
  MatrixToAffineParameters(m, tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz);
  EXPECT_DOUBLE_EQ( 0.0, tx);
  EXPECT_DOUBLE_EQ( 0.0, ty);
  EXPECT_DOUBLE_EQ( 0.0, tz);
  EXPECT_DOUBLE_EQ( 0.7, rx);
  EXPECT_DOUBLE_EQ(-1.1, ry);
  EXPECT_DOUBLE_EQ( 3.0, rz);
  EXPECT_DOUBLE_EQ( 1.0, sx);
  EXPECT_DOUBLE_EQ( 1.0, sy);
  EXPECT_DOUBLE_EQ( 1.0, sz);
  EXPECT_DOUBLE_EQ( 0.0, sxy);
  EXPECT_DOUBLE_EQ( 0.0, sxz);
  EXPECT_DOUBLE_EQ( 0.0, syz);
}

// -----------------------------------------------------------------------------
TEST(Matrix, ScalingToAffineParameters)
{
  Matrix m;
  double tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz;
  AffineParametersToMatrix(.0, .0, .0, .0, .0, .0, 4.0, 1.7, 0.3, .0, .0, .0, m);
  MatrixToAffineParameters(m, tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz);
  EXPECT_DOUBLE_EQ( 0.0, tx);
  EXPECT_DOUBLE_EQ( 0.0, ty);
  EXPECT_DOUBLE_EQ( 0.0, tz);
  EXPECT_DOUBLE_EQ( 0.0, rx);
  EXPECT_DOUBLE_EQ( 0.0, ry);
  EXPECT_DOUBLE_EQ( 0.0, rz);
  EXPECT_DOUBLE_EQ( 4.0, sx);
  EXPECT_DOUBLE_EQ( 1.7, sy);
  EXPECT_DOUBLE_EQ( 0.3, sz);
  EXPECT_DOUBLE_EQ( 0.0, sxy);
  EXPECT_DOUBLE_EQ( 0.0, sxz);
  EXPECT_DOUBLE_EQ( 0.0, syz);
}

// -----------------------------------------------------------------------------
TEST(Matrix, ShearingToAffineParameters)
{
  Matrix m;
  double tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz;
  AffineParametersToMatrix(.0, .0, .0, .0, .0, .0, 1.0, 1.0, 1.0, .7, -1.1, 1.5, m);
  MatrixToAffineParameters(m, tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz);
  EXPECT_DOUBLE_EQ( 0.0, tx);
  EXPECT_DOUBLE_EQ( 0.0, ty);
  EXPECT_DOUBLE_EQ( 0.0, tz);
  EXPECT_DOUBLE_EQ( 0.0, rx);
  EXPECT_DOUBLE_EQ( 0.0, ry);
  EXPECT_DOUBLE_EQ( 0.0, rz);
  EXPECT_DOUBLE_EQ( 1.0, sx);
  EXPECT_DOUBLE_EQ( 1.0, sy);
  EXPECT_DOUBLE_EQ( 1.0, sz);
  EXPECT_DOUBLE_EQ( 0.7, sxy);
  EXPECT_DOUBLE_EQ(-1.1, sxz);
  EXPECT_DOUBLE_EQ( 1.5, syz);
}

// -----------------------------------------------------------------------------
TEST(Matrix, TranslationAndRotationToAffineParameters)
{
  Matrix m;
  double tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz;
  AffineParametersToMatrix(-1.0, 0.5, 4.7, .7, -1.1, 3.0, 1.0, 1.0, 1.0, .0, .0, .0, m);
  MatrixToAffineParameters(m, tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz);
  EXPECT_DOUBLE_EQ(-1.0, tx);
  EXPECT_DOUBLE_EQ( 0.5, ty);
  EXPECT_DOUBLE_EQ( 4.7, tz);
  EXPECT_DOUBLE_EQ( 0.7, rx);
  EXPECT_DOUBLE_EQ(-1.1, ry);
  EXPECT_DOUBLE_EQ( 3.0, rz);
  EXPECT_DOUBLE_EQ( 1.0, sx);
  EXPECT_DOUBLE_EQ( 1.0, sy);
  EXPECT_DOUBLE_EQ( 1.0, sz);
  EXPECT_DOUBLE_EQ( 0.0, sxy);
  EXPECT_DOUBLE_EQ( 0.0, sxz);
  EXPECT_DOUBLE_EQ( 0.0, syz);
}

// -----------------------------------------------------------------------------
TEST(Matrix, ScalingAndShearingToAffineParameters)
{
  Matrix m;
  double tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz;
  AffineParametersToMatrix(.0, .0, .0, .0, .0, .0, 4.0, 1.7, 0.3, .7, -1.1, 1.5, m);
  MatrixToAffineParameters(m, tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz);
  EXPECT_DOUBLE_EQ( 0.0, tx);
  EXPECT_DOUBLE_EQ( 0.0, ty);
  EXPECT_DOUBLE_EQ( 0.0, tz);
  EXPECT_DOUBLE_EQ( 0.0, rx);
  EXPECT_DOUBLE_EQ( 0.0, ry);
  EXPECT_DOUBLE_EQ( 0.0, rz);
  EXPECT_DOUBLE_EQ( 4.0, sx);
  EXPECT_DOUBLE_EQ( 1.7, sy);
  EXPECT_DOUBLE_EQ( 0.3, sz);
  EXPECT_DOUBLE_EQ( 0.7, sxy);
  EXPECT_DOUBLE_EQ(-1.1, sxz);
  EXPECT_DOUBLE_EQ( 1.5, syz);
}

// -----------------------------------------------------------------------------
TEST(Matrix, TranslationRotationScalingAndShearingToAffineParameters)
{
  Matrix m;
  double tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz;
  AffineParametersToMatrix(-1.0, 0.5, 4.7, .7, -1.1, 3.0, 4.0, 1.7, 0.3, .7, -1.1, 1.5, m);
  MatrixToAffineParameters(m, tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz);
  EXPECT_DOUBLE_EQ(-1.0, tx);
  EXPECT_DOUBLE_EQ( 0.5, ty);
  EXPECT_DOUBLE_EQ( 4.7, tz);
  EXPECT_DOUBLE_EQ( 0.7, rx);
  EXPECT_DOUBLE_EQ(-1.1, ry);
  EXPECT_DOUBLE_EQ( 3.0, rz);
  EXPECT_DOUBLE_EQ( 4.0, sx);
  EXPECT_DOUBLE_EQ( 1.7, sy);
  EXPECT_DOUBLE_EQ( 0.3, sz);
  EXPECT_DOUBLE_EQ( 0.7, sxy);
  EXPECT_DOUBLE_EQ(-1.1, sxz);
  EXPECT_DOUBLE_EQ( 1.5, syz);
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
TEST(Matrix, IO)
{
  // Using STL streams
  {
    // Write matrix to binary file
    ofstream ofs("testMatrix1.bin", ios::binary);
    ASSERT_TRUE(ofs.is_open());
    Matrix m1 = AffineTestMatrix();
    ofs << m1;
    ASSERT_FALSE(ofs.fail());
    ofs.close();
    ASSERT_FALSE(ofs.is_open());
    // Read binary file and check first few binary bytes to see if data
    // is in big-endian order as expected
    ifstream tmp("testMatrix1.bin", ios::binary);
    ASSERT_TRUE(tmp.is_open());
    char buffer[32];
    tmp.read(buffer, 32);
    tmp.close();
    ASSERT_TRUE(strncmp(buffer, "irtkMatrix 4 x 4\n", 17) == 0);
    EXPECT_EQ(buffer[17], (char) 63); // 3F
    EXPECT_EQ(buffer[18], (char)234); // EA
    EXPECT_EQ(buffer[19], (char)109); // 6D
    // Read matrix from binary file stream
    ifstream ifs("testMatrix1.bin", ios::binary);
    ASSERT_TRUE(ifs.is_open());
    Matrix m2;
    ifs >> m2;
    ASSERT_FALSE(ifs.fail());
    ifs.close();
    // Check that read matrix is the same that we saved before
    ASSERT_TRUE(m1 == m2);
  }
  // Using MIRTK streams
  {
    // Write matrix to binary file
    Cofstream ofs("testMatrix2.bin");
    Matrix m1 = AffineTestMatrix();
    ofs << m1;
    ofs.Close();
    // Read binary file and check first few binary bytes to see if data
    ifstream tmp("testMatrix2.bin", ios::binary);
    ASSERT_TRUE(tmp.is_open());
    char buffer[32];
    tmp.read(buffer, 32);
    tmp.close();
    ASSERT_TRUE(strncmp(buffer, "irtkMatrix", 11) == 0);
    EXPECT_EQ(buffer[14], (char)4);
    EXPECT_EQ(buffer[18], (char)4);
    EXPECT_EQ(buffer[19], (char) 63); // 3F
    EXPECT_EQ(buffer[20], (char)234); // EA
    EXPECT_EQ(buffer[21], (char)109); // 6D
    // Read matrix from binary file again
    Cifstream ifs("testMatrix2.bin");
    Matrix m2;
    ifs >> m2;
    ifs.Close();
    ASSERT_TRUE(m1 == m2);
  }
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
