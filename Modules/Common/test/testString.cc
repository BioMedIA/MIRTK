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

#include "gtest/gtest.h"

#include "mirtk/String.h"
using namespace mirtk;


// =============================================================================
// String conversion
// =============================================================================

// -----------------------------------------------------------------------------
TEST(String, FromString)
{
  {
    bool b;
    EXPECT_FALSE(FromString("1.0", b));
    EXPECT_FALSE(FromString("a", b));
    EXPECT_TRUE(FromString("1", b));
    EXPECT_TRUE(b);
    EXPECT_TRUE(FromString("0", b));
    EXPECT_FALSE(b);
    EXPECT_TRUE(FromString("On", b));
    EXPECT_TRUE(b);
  }
  {
    double v;
    EXPECT_TRUE(FromString("1.0", v));
    EXPECT_EQ(1.0, v);
  }
}

// =============================================================================
// String processing
// =============================================================================

// -----------------------------------------------------------------------------
TEST(String, ToLower)
{
  EXPECT_EQ(string(""),      ToLower(string("")));
  EXPECT_EQ(string("test"),  ToLower(string("test")));
  EXPECT_EQ(string("test"),  ToLower(string("Test")));
  EXPECT_EQ(string("test"),  ToLower(string("TeSt")));
  EXPECT_EQ(string("test"),  ToLower(string("TEST")));
  EXPECT_EQ(string("test"),  ToLower(string("tesT")));
  EXPECT_EQ(string(" test"), ToLower(string(" tesT")));
}

// -----------------------------------------------------------------------------
TEST(String, ToUpper)
{
  EXPECT_EQ(string(""),      ToUpper(string("")));
  EXPECT_EQ(string("TEST"),  ToUpper(string("TEST")));
  EXPECT_EQ(string("TEST"),  ToUpper(string("tEST")));
  EXPECT_EQ(string("TEST"),  ToUpper(string("tEsT")));
  EXPECT_EQ(string("TEST"),  ToUpper(string("test")));
  EXPECT_EQ(string("TEST"),  ToUpper(string("TESt")));
  EXPECT_EQ(string(" TEST"), ToUpper(string(" TESt")));
}

// -----------------------------------------------------------------------------
TEST(String, Trim)
{
  EXPECT_EQ(string("test"), Trim("   \ttest  \t  "));
  EXPECT_EQ(string("test"), Trim("test  \t  \n"));
  EXPECT_EQ(string("test"), Trim(" \n  \ttest"));
}

// -----------------------------------------------------------------------------
TEST(String, Split)
{
  const int n = 5;
  const char *delim[n] = {" ", ",", "\t", ";", "_@#$"};
  {
    Array<string> parts = Split("a b c d", " ");
    ASSERT_EQ(4u, parts.size());
    EXPECT_EQ(string("a"), parts[0]);
    EXPECT_EQ(string("b"), parts[1]);
    EXPECT_EQ(string("c"), parts[2]);
    EXPECT_EQ(string("d"), parts[3]);
  }
  {
    for (int i = 0; i < n; ++i) {
      string str = string("a") + delim[i] + delim[i] + "b" + delim[i] + "c";
      for (int j = 0; j < 5; ++j) {
        str += delim[i];
      }
      str += "d";
      Array<string> parts = Split(str, delim[i]);
      ASSERT_EQ(9u, parts.size());
      EXPECT_EQ(string("a"), parts[0]);
      EXPECT_EQ(string(""),  parts[1]);
      EXPECT_EQ(string("b"), parts[2]);
      EXPECT_EQ(string("c"), parts[3]);
      EXPECT_EQ(string(""),  parts[4]);
      EXPECT_EQ(string(""),  parts[5]);
      EXPECT_EQ(string(""),  parts[6]);
      EXPECT_EQ(string(""),  parts[7]);
      EXPECT_EQ(string("d"), parts[8]);
    }
  }
  {
    const bool discard_empty = true;
    for (int i = 0; i < n; ++i) {
      string str = string("a") + delim[i] + delim[i] + "b" + delim[i] + "c";
      for (int j = 0; j < 5; ++j) {
        str += delim[i];
      }
      str += "d";
      Array<string> parts = Split(str, delim[i], 0, discard_empty);
      ASSERT_EQ(4u, parts.size());
      EXPECT_EQ(string("a"), parts[0]);
      EXPECT_EQ(string("b"), parts[1]);
      EXPECT_EQ(string("c"), parts[2]);
      EXPECT_EQ(string("d"), parts[3]);
    }
  }
  {
    const bool discard_empty = true;
    for (int i = 0; i < n; ++i) {
      string str = string("a") + delim[i] + delim[i] + "b" + delim[i] + "c";
      for (int j = 0; j < 5; ++j) {
        str += delim[i];
      }
      str += "d";
      Array<string> parts = Split(str, delim[i], -1, discard_empty);
      ASSERT_EQ(1u, parts.size());
      EXPECT_EQ(string("d"), parts[0]);
    }
  }
  {
    const bool discard_empty = true;
    for (int i = 0; i < n; ++i) {
      string str = string("a") + delim[i] + delim[i] + "b" + delim[i] + "c";
      for (int j = 0; j < 5; ++j) {
        str += delim[i];
      }
      str += "d";
      Array<string> parts = Split(str, delim[i], -4, discard_empty);
      ASSERT_EQ(4u, parts.size());
      EXPECT_EQ(string("a"), parts[0]);
      EXPECT_EQ(string("b"), parts[1]);
      EXPECT_EQ(string("c"), parts[2]);
      EXPECT_EQ(string("d"), parts[3]);
    }
  }
  {
    const bool discard_empty = false;
    const bool handle_quotes = true;
    auto parts = Split(" \t 'a'  \"b\"b  \tc \t\"d\"", " ", 0, discard_empty, handle_quotes);
    ASSERT_EQ(5u, parts.size());
    EXPECT_TRUE(parts[0].empty());
    EXPECT_EQ(string("\t"), parts[1]);
    EXPECT_EQ(string("a"), parts[2]);
    EXPECT_TRUE(parts[3].empty());
    EXPECT_EQ(string("b\"b  \tc \t\"d"), parts[4]);
  }
  {
    const bool discard_empty = false;
    const bool handle_quotes = true;
    auto parts = Split(" \t 'a'  \"b\"b  \tc \t\"d\"", " ", -2, discard_empty, handle_quotes);
    EXPECT_TRUE(parts[0].empty());
    EXPECT_EQ(string("b\"b  \tc \t\"d"), parts[1]);
  }
  {
    const bool discard_empty = true;
    const bool handle_quotes = true;
    auto parts = Split(" \t 'a'  \"b\"b  \tc \t\"d\"", " ", -2, discard_empty, handle_quotes);
    EXPECT_EQ(string("a"), parts[0]);
    EXPECT_EQ(string("b\"b  \tc \t\"d"), parts[1]);
  }
  {
    const bool discard_empty = true;
    const bool handle_quotes = true;
    auto parts = Split(" first 'sec ond '  t \t\"fo\"ur  th\"", " ", 0, discard_empty, handle_quotes);
    ASSERT_EQ(4u, parts.size());
    EXPECT_EQ(string("first"), parts[0]);
    EXPECT_EQ(string("sec ond "), parts[1]);
    EXPECT_EQ(string("t"), parts[2]);
    EXPECT_EQ(string("fo\"ur  th"), parts[3]);
  }
  {
    const bool discard_empty = true;
    const bool handle_quotes = true;
    auto parts = Split(" first,'sec o,nd ' , t,\t \"fo\"ur  th\"\n", ",", 0, discard_empty, handle_quotes);
    ASSERT_EQ(4u, parts.size());
    EXPECT_EQ(string(" first"), parts[0]);
    EXPECT_EQ(string("sec o,nd "), parts[1]);
    EXPECT_EQ(string(" t"), parts[2]);
    EXPECT_EQ(string("fo\"ur  th"), parts[3]);
  }
}

// -----------------------------------------------------------------------------
TEST(String, StandardUnits)
{
  EXPECT_EQ(string("vox"), StandardUnits("voxels"));
  EXPECT_EQ(string("vox"), StandardUnits("VOX"));
  EXPECT_EQ(string("vox"), StandardUnits("VOXEL"));
  EXPECT_EQ(string("%"),   StandardUnits("percentage"));
  EXPECT_EQ(string("mm"),  StandardUnits("millimeters"));
  EXPECT_EQ(string("mm"),  StandardUnits("MM"));
  EXPECT_EQ(string("rel"), StandardUnits("relative"));
  EXPECT_EQ(string("abs"), StandardUnits("ABSOLUTE"));
  EXPECT_EQ(string("abs"), StandardUnits("absolute"));
}

// -----------------------------------------------------------------------------
TEST(String, ParameterUnits)
{
  {
    string name;
    EXPECT_EQ(string("mm"), ParameterUnits("Resolution[mm]  \n", &name));
    EXPECT_EQ("Resolution", name);
  }
  EXPECT_EQ(string("mm"), ParameterUnits("Resolution [mm]"));
  EXPECT_EQ(string("mm"), ParameterUnits("Resolution [MM]"));
  {
    string name;
    EXPECT_EQ(string("mm"), ParameterUnits("  Resolution    [MM]", &name));
    EXPECT_EQ("Resolution", name);
  }
  {
    string name;
    EXPECT_EQ(string(""), ParameterUnits("Resolution", &name));
    EXPECT_EQ("Resolution", name);
  }
  {
    string name;
    EXPECT_EQ(string(""), ParameterUnits("Resolution of image 1", &name));
    EXPECT_EQ("Resolution of image 1", name);
  }
  EXPECT_EQ(string("rel"), ParameterUnits("Resolution [relative]  \t"));
}

// -----------------------------------------------------------------------------
TEST(String, ValueUnits)
{
  EXPECT_EQ(string("mm"), ValueUnits("1 mm"));
  EXPECT_EQ(string("mm"), ValueUnits("1mm"));
  EXPECT_EQ(string("mm"), ValueUnits("1  [mm]"));
  {
    string value;
    EXPECT_EQ(string("mm"), ValueUnits("  1[mm]    \n", &value));
    EXPECT_EQ(string("1"), value);
  }
  {
    string value;
    EXPECT_EQ(string("mm"), ValueUnits("  1 2 3mm", &value));
    EXPECT_EQ(string("1 2 3"), value);
  }
  {
    string value;
    EXPECT_EQ(string("mm"), ValueUnits("  1.5mm", &value));
    EXPECT_EQ(string("1.5"), value);
  }
  {
    string value;
    EXPECT_EQ(string("mm"), ValueUnits("  1 2 3  [mm]", &value));
    EXPECT_EQ(string("1 2 3"), value);
  }
  EXPECT_EQ(string(""), ValueUnits("1"));
  EXPECT_EQ(string(""), ValueUnits("foo mm"));
  EXPECT_EQ(string("mm"), ValueUnits("foo [mm]\n"));
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
