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

#include "mirtk/Polyhedron.h"

#include "mirtk/Assert.h"
#include "mirtk/Parallel.h"

#include "vtkPolyData.h"
#include "vtkCellArray.h"


namespace mirtk {


// =============================================================================
// Auxiliaries
// =============================================================================

namespace PolyhedronUtils {


// -----------------------------------------------------------------------------
/// Compute volume of polyhedron
struct ComputeVolume
{
  vtkPolyData *_DataSet;
  double       _Sum;

  ComputeVolume() : _Sum(.0) {}

  ComputeVolume(const ComputeVolume &other, split)
  :
    _DataSet(other._DataSet), _Sum(.0)
  {}

  void join(const ComputeVolume &other)
  {
    _Sum += other._Sum;
  }

  void operator ()(const blocked_range<vtkIdType> &re)
  {
    double v1[3], v2[3], v3[3], det, height;
    vtkIdType npts, *pts;
    vtkCellArray * const polys = _DataSet->GetPolys();
    for (vtkIdType cellId = re.begin(); cellId != re.end(); ++cellId) {
      // Get position of triangle vertices
      polys->GetCell(cellId, npts, pts);
      if (npts == 0) continue;
      mirtkAssert(npts == 3, "mesh is triangulated");
      _DataSet->GetPoint(pts[0], v1);
      _DataSet->GetPoint(pts[1], v2);
      _DataSet->GetPoint(pts[2], v3);
      // Twice the area of the projection onto the x-y plane
      det = ((v2[1] - v3[1]) * (v1[0] - v3[0]) -
             (v2[0] - v3[0]) * (v1[1] - v3[1]));
      // Three times the average height
      height = v1[2] + v2[2] + v3[2];
      _Sum += det * height;
    }
  }

  double Value() const
  {
    return _Sum / 6.0;
  }
};

// -----------------------------------------------------------------------------
/// Compute winding number of polyhedron around a given point
/*
Robust point-in-polyhedron test description and implementation by Mark Dickinson.
https://github.com/mdickinson/polyhedron/blob/cd7361bcee8cbd9ef2e0aac58a7ce59bf9a52c4f/polyhedron.py

Given an closed, oriented surface in R^3 described by a triangular mesh, the
code below gives a robust algorithm for determining whether a given point is
inside, on the boundary of, or outside, the surface.  The algorithm should give
correct results even in degenerate cases, and applies to disconnected
polyhedra, non simply-connected surfaces, and so on.  There are no requirements
for the surface to be convex, simple, connected or simply-connected.

More precisely, we give a method for computing the *winding number* of a closed
oriented surface S around a point O that doesn't lie on S.  Roughly speaking,
the winding number of the closed oriented surface S around a point O not on S
is the number of times that the surface encloses that point; for a simple
outward-oriented surface (like that of a convex polyhedron, for example), the
winding number will be 1 for points inside the surface and 0 for points
outside.

For a precise definition of winding number, we can turn to algebraic topology:
our oriented surface is presented as a collection of combinatorial data
defining abstract vertices, edges and triangles, together with a mapping of
those vertices to R^3.  The combinatorial data describe a simplicial complex C,
and assuming that O doesn't lie on the surface, the mapping of the vertices to
R^3 gives a continuous map from the geometric realization of C to R^3 - {O}.
This in turn induces a map on second homology groups:

   H^2(C, Z) -> H^2(R^3 - {O}, Z)

and by taking the usual right-handed orientation in R^3 we identify H^2(R^3 -
{O}, Z) with Z.  The image of [S] under this map gives the winding number.  In
particular, the well-definedness of the winding number does not depend on
topological properties of the embedding: it doesn't matter if the surface is
self-intersecting, or has degenerate triangles.  The only condition is that O
does not lie on the surface S.

Algorithm
---------

The algorithm is based around the usual method of ray-casting: we take a
vertical line L through O and count the intersections of this line with the
triangles of the surface, keeping track of orientations as we go.  Let's ignore
corner cases for a moment and assume that:

(1) O does not lie on the surface, and
(2) for each triangle T (thought of as a closed subset of R^3) touched by
    our vertical line L, L meets the interior of T in exactly one point Q

Then there are four possibilities for each such triangle T:
1. T lies *above* O and is oriented *upwards* (*away* from O).
2. T lies *above* O and is oriented *downwards* (*towards* O).
3. T lies *below* O and is oriented *downwards* (*away* from O).
4. T lies *below* O and is oriented *upwards* (*towards* O).

Let's write N1, N2, N3 and N4 for the counts of triangles satisfying conditions
1, 2, 3 and 4 respectively.  Since we have a closed surface, these numbers
are not independent; they satisfy the relation:

    N1 + N4 == N2 + N3

That is, the number of upward-facing triangles must match the number of
downward-facing triangles.  The winding number w is then given by:

    w = N1 - N2 == N3 - N4

In the code below, we simply compute 2*w = (N1 + N3) - (N2 + N4), so each
triangle oriented away from O contributes 1 to 2w, while each triangle oriented
towards O contributes -1.

Making the algorithm robust
---------------------------

Now we describe how to augment the basic algorithm described above to include:

- correct treatment of corner cases (vertical triangles, cases where L meets an
  edge or vertex directly, etc.)
- detection of the case where the point lies directly on the surface.

It turns out that to make the algorithm robust, all we need to do is be careful
and consistent about classifying vertices, edges and triangles.  We do this as
follows:

- Each vertex of the surface that's not equal to O is considered *positive* if
  its coordinates are lexicographically greater than O, and *negative*
  otherwise.
- For an edge PQ of the surface that's not collinear with O, we first describe
  the classification in the case that P is negative and Q is positive, and
  then extend to arbitrary PQ.
  For P negative and Q positive, there are two cases:
  1. P and Q have distinct x coordinates.  In that case we classify the edge
     PQ by its intersection with the plane passing through O and parallel
     to the yz-plane: the edge is *positive* if the intersection point is
     positive, and *negative* otherwise.
  2. P and Q have the same x coordinate, in which case they must have
     distinct y coordinates.  (If the x and the y coordinates both match
     then PQ passes through O.)  We classify by the intersection of PQ
     with the line parallel to the y-axis through O.
  For P positive and Q negative, we classify as above but reverse the sign.
  For like-signed P and Q, the classification isn't used.
  Computationally, in case 1 above, the y-coordinate of the intersection
  point is:

      Py + (Qy - Py) * (Ox - Px) / (Qx - Px)

  and this is greater than Oy iff

      (Py - Oy) * (Qx - Ox) - (Px - Ox) * (Qy - Oy)

  is positive, so the sign of the edge is the sign of the above expression.
  Similarly, if this quantity is zero then we need to look at the z-coordinate
  of the intersection, and the sign of the edge is given by

      (Pz - Oz) * (Qx - Ox) - (Px - Ox) * (Qz - Oz)

  In case 2, both of the above quantities are zero, and the sign of the edge is
  the sign of

      (Pz - Oz) * (Qy - Oy) - (Py - Oy) * (Qz - Oz)

  Another way to look at this: if P, Q and O are not collinear then the
  matrix

   ( Px Qx Ox )
   ( Py Qy Ox )
   ( Pz Qz Ox )
   (  1  1  1 )

  has rank 3.  It follows that at least one of the three 3x3 minors

   | Px Qx Ox |  | Px Qx Ox |  | Py Qy Oy |
   | Py Qy Oy |  | Pz Qz Oz |  | Pz Qz Oz |
   |  1  1  1 |  |  1  1  1 |  |  1  1  1 |

  is nonzero.  We define the sign of PQ to be the *negative* of the sign of the
  first nonzero minor in that list.

- Each triangle PQR of the surface that's not coplanar with O is considered
  *positive* if its normal points away from O, and *negative* if its normal
  points towards O.
  Computationally, the sign of the triangle PQR is the sign of the 4x4
  determinant

    | Px Qx Rx Ox |
    | Py Qy Ry Oy |
    | Pz Qz Rz Oz |
    |  1  1  1  1 |

  or equivalently of the 3x3 determinant

    | Px-Ox Qx-Ox Rx-Ox |
    | Py-Oy Qy-Oy Ry-Oy |
    | Pz-Oz Qz-Oz Rz-Oz |

Now to compute the contribution of any given triangle to the total winding
number:

1. Classify the vertices of the triangle.  At the same time, we can check that
   none of the vertices is equal to O.  If all vertices have the same sign,
   then the winding number contribution is zero.
2. Assuming that the vertices do not all have the same sign, two of the three
   edges connect two differently-signed vertices.  Classify both those edges
   (and simultaneously check that they don't pass through O).  If the edges
   have opposite classification, then the winding number contribution is zero.
3. Now two of the edges have the same sign: classify the triangle itself.  If
   the triangle is positive it contributes 1/2 to the winding number total; if
   negative it contributes -1/2.  In practice we count contributions of 1 and
   -1, and halve the total at the end.

Note that an edge between two like-signed vertices can never pass through O, so
there's no need to check the third edge in step 2.  Similarly, a triangle whose
edge-cycle is trivial can't contain O in its interior.

To understand what's going on above, it's helpful to step into the world of
homology again. The homology of R^3 - {O} can be identified with that of the
two-sphere S^2 by deformation retract, and we can decompose the two-sphere as a
CW complex consisting of six cells, as follows:

* 0-cells B and F, where B = (-1, 0, 0) and F = (1, 0, 0)
* 1-cells L and R, where
     L = {(cos t, sin t, 0) | -pi <= t <= 0 }
     R = {(cos t, sin t, 0) | 0 <= t <= pi }
* 2-cells U and D, where U is the top half of the sphere (z >= 0)
  and D is the bottom half (z <= 0), both oriented outwards.

And the homology of the CW complex is now representable in terms of cellular
homology:
               d               d
  Z[U] + Z[D] --> Z[L] + Z[R] --> Z[B] + Z[F]

with boundary maps given by:

  d[U] = [L] + [R]; d[D] = -[L] - [R]
  d[R] = [B] - [F]; d[L] = [F] - [B]

Now the original map C -> R^3 - {O} from the geometric realization of the
simplicial complex is homotopic to a map C -> S^2 that sends:

* each positive vertex to F and each negative vertex to B
* each edge with boundary [F] - [B] to L if the edge is negative, and -R if the
  edge is positive
* each edge with boundary [B] - [F] to R if the edge is positive, and -L if the
  edge is negative
* all other edges to 0
* each triangle whose boundary is [L] + [R] to either U or -D,
  depending on whether the triangle is positive or negative
* each triangle whose boundary is -[L] - [R] to either D or -U,
  depending on whether the triangle is positive or negative
* all other triangles to 0

Mapping all of the triangles in the surface this way, and summing the results
in second homology, we end up with (winding number)*([U] + [D]).
*/
struct ComputeWindingNumber
{
  vtkPolyData *_DataSet;
  double       _X, _Y, _Z;
  int          _Sum;

  /// Default constructor
  ComputeWindingNumber() : _Sum(0) {}

  /// Split constructor
  ComputeWindingNumber(const ComputeWindingNumber &other, split)
  :
    _DataSet(other._DataSet), _X(other._X), _Y(other._Y), _Z(other._Z), _Sum(0)
  {}

  // ---------------------------------------------------------------------------
  /// Join contributions of threads
  void join(const ComputeWindingNumber &other)
  {
    _Sum += other._Sum;
  }

  // ---------------------------------------------------------------------------
  /// Return 1 if x is positive, -1 if it's negative, and 0 if it's zero.
  inline int sign(double x) const
  {
    return static_cast<int>(x > 0) - static_cast<int>(x < 0);
  }

  // ---------------------------------------------------------------------------
  /// Sign of the vertex P with respect to O, as defined above.
  inline int vertex_sign(double P[3]) const
  {
    int              result = sign(P[0] - _X);
    if (result == 0) result = sign(P[1] - _Y);
    if (result == 0) result = sign(P[2] - _Z);
    if (result == 0) {
      //throw new invalid_argument("vertex coincides with origin");
    }
    return result;
  }

  // ---------------------------------------------------------------------------
  /// Sign of the edge PQ with respect to O, as defined above.
  inline int edge_sign(double P[3], double Q[3]) const
  {
    int result = sign((P[1] - _Y) * (Q[0] - _X) - (P[0] - _X) * (Q[1] - _Y));
    if (result == 0) {
      result = sign((P[2] - _Z) * (Q[0] - _X) - (P[0] - _X) * (Q[2] - _Z));
      if (result == 0) {
        result = sign((P[2] - _Z) * (Q[1] - _Y) - (P[1] - _Y) * (Q[2] - _Z));
        if (result == 0) {
          //throw new invalid_argument("vertices collinear with origin");
        }
      }
    }
    return result;
  }

  // ---------------------------------------------------------------------------
  /// Sign of the triangle PQR with respect to O, as defined above.
  inline int triangle_sign(double P[3], double Q[3], double R[3]) const
  {
    const double m1_0 = P[0] - _X;
    const double m1_1 = P[1] - _Y;
    const double m2_0 = Q[0] - _X;
    const double m2_1 = Q[1] - _Y;
    const double m3_0 = R[0] - _X;
    const double m3_1 = R[1] - _Y;
    int result = sign((m1_0 * m2_1 - m1_1 * m2_0) * (R[2] - _Z) +
                      (m2_0 * m3_1 - m2_1 * m3_0) * (P[2] - _Z) +
                      (m3_0 * m1_1 - m3_1 * m1_0) * (Q[2] - _Z));
    if (result == 0) {
      //throw new invalid_argument("vertices coplanar with origin");
    }
    return result;
  }

  // ---------------------------------------------------------------------------
  /// Return the contribution of this triangle to the winding number.
  /// \throws invalid_argument if the face contains the origin.
  inline int triangle_chain(double v1[3], double v2[3], double v3[3]) const
  {
    const int v1sign = vertex_sign(v1);
    const int v2sign = vertex_sign(v2);
    const int v3sign = vertex_sign(v3);

    int face_boundary = 0;
    if (v1sign != v2sign) face_boundary += edge_sign(v1, v2);
    if (v2sign != v3sign) face_boundary += edge_sign(v2, v3);
    if (v3sign != v1sign) face_boundary += edge_sign(v3, v1);
    if (face_boundary == 0) return 0;

    return triangle_sign(v1, v2, v3);
  }

  // ---------------------------------------------------------------------------
  /// Add contributions of triangles in range to winding number
  void operator ()(const blocked_range<vtkIdType> &re)
  {
    double v1[3], v2[3], v3[3];
    vtkIdType npts, *pts;
    vtkCellArray * const polys = _DataSet->GetPolys();
    for (vtkIdType cellId = re.begin(); cellId != re.end(); ++cellId) {
      polys->GetCell(cellId, npts, pts);
      if (npts == 0) continue;
      mirtkAssert(npts == 3, "mesh is triangulated");
      _DataSet->GetPoint(pts[0], v1);
      _DataSet->GetPoint(pts[1], v2);
      _DataSet->GetPoint(pts[2], v3);
      _Sum += triangle_chain(v1, v2, v3);
    }
  }

  // ---------------------------------------------------------------------------
  /// Get winding number
  int Value() const
  {
    return _Sum / 2;
  }
};


} // namespace PolyhedronUtils

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
Polyhedron::Polyhedron(vtkPolyData *polydata)
:
  _DataSet(polydata)
{
}

// -----------------------------------------------------------------------------
Polyhedron::Polyhedron(const Polyhedron &other)
:
  Object(other),
  _DataSet(other._DataSet)
{
}

// -----------------------------------------------------------------------------
Polyhedron &Polyhedron::operator =(const Polyhedron &other)
{
  if (this != &other) {
    Object::operator =(other);
    _DataSet = other._DataSet;
  }
  return *this;
}

// -----------------------------------------------------------------------------
Polyhedron::~Polyhedron()
{
}

// =============================================================================
// Properties
// =============================================================================

// -----------------------------------------------------------------------------
double Polyhedron::Volume(vtkPolyData *polydata)
{
  PolyhedronUtils::ComputeVolume vol;
  vol._DataSet = polydata;
  blocked_range<vtkIdType> cellIds(0, polydata->GetPolys()->GetNumberOfCells());
  parallel_reduce(cellIds, vol);
  return vol.Value();
}

// -----------------------------------------------------------------------------
int Polyhedron::WindingNumber(vtkPolyData *polydata, double x, double y, double z)
{
  PolyhedronUtils::ComputeWindingNumber wn;
  wn._DataSet = polydata;
  wn._X = x, wn._Y = y, wn._Z = z;
  blocked_range<vtkIdType> cellIds(0, polydata->GetPolys()->GetNumberOfCells());
  parallel_reduce(cellIds, wn);
  return wn.Value();
}


} // namespace mirtk
