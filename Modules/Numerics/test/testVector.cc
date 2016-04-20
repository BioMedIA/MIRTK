/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#include "mirtk/Vector.h"
using namespace mirtk;


// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
Vector IncreasingTestVector(int n = 10)
{
  Vector v(n);
  for (int i = 0; i < n; ++i) v(i) = double(i+1);
  return v;
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
TEST(Vector, IO)
{
  // Using STL streams
  {
    // Write vector to binary file
    ofstream ofs("testVector1.bin", ios::binary);
    ASSERT_TRUE(ofs.is_open());
    Vector v1 = IncreasingTestVector(10);
    ofs << v1;
    ASSERT_FALSE(ofs.fail());
    ofs.close();
    ASSERT_FALSE(ofs.is_open());
    // Read binary file and check first few binary bytes to see if data
    // is in big-endian order as expected
    ifstream tmp("testVector1.bin", ios::binary);
    ASSERT_TRUE(tmp.is_open());
    char buffer[64];
    tmp.read(buffer, 64);
    tmp.close();
    ASSERT_TRUE(strncmp(buffer, "irtkVector 10\n", 14) == 0);
    EXPECT_EQ((char)63,  buffer[14]);
    EXPECT_EQ((char)240, buffer[15]);
    EXPECT_EQ((char)64,  buffer[22]);
    EXPECT_EQ((char)0,   buffer[23]);
    EXPECT_EQ((char)64,  buffer[30]);
    EXPECT_EQ((char)8,   buffer[31]);
    // Read vector from binary file stream
    ifstream ifs("testVector1.bin", ios::binary);
    ASSERT_TRUE(ifs.is_open());
    Vector v2;
    ifs >> v2;
    ASSERT_FALSE(ifs.fail());
    ifs.close();
    // Check that read vector is the same that we saved before
    ASSERT_TRUE(v1 == v2);
  }
  // Using MIRTK streams
  {
    // Write vector to binary file
    Cofstream ofs("testVector2.bin");
    Vector v1 = IncreasingTestVector(10);
    ofs << v1;
    ofs.Close();
    // Read binary file and check first few binary bytes to see if data
    ifstream tmp("testVector2.bin", ios::binary);
    ASSERT_TRUE(tmp.is_open());
    char buffer[64];
    tmp.read(buffer, 64);
    tmp.close();
    ASSERT_TRUE(strncmp(buffer, "irtkVector", 11) == 0);
    EXPECT_EQ((char)0,   buffer[11]);
    EXPECT_EQ((char)0,   buffer[12]);
    EXPECT_EQ((char)0,   buffer[13]);
    EXPECT_EQ((char)10,  buffer[14]);
    EXPECT_EQ((char)63,  buffer[15]);
    EXPECT_EQ((char)240, buffer[16]);
    EXPECT_EQ((char)64,  buffer[23]);
    EXPECT_EQ((char)0,   buffer[24]);
    EXPECT_EQ((char)64,  buffer[31]);
    EXPECT_EQ((char)8,   buffer[32]);
    // Read vector from binary file again
    Cifstream ifs("testVector2.bin");
    Vector v2;
    ifs >> v2;
    ifs.Close();
    ASSERT_TRUE(v1 == v2);
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
