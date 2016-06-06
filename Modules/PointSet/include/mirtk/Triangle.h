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

#ifndef MIRTK_Triangle_H
#define MIRTK_Triangle_H

#include "mirtk/Memory.h"
#include "mirtk/Math.h"
#include "mirtk/VtkMath.h"

#include "vtkLine.h"
#include "vtkPlane.h"
#include "vtkTriangle.h"


namespace mirtk {


/**
 * Auxiliary class/static helper functions for triangulated surface meshes
 */
class Triangle
{
public:

  /// Compute normal direction of triangle
  ///
  /// \param[in]  a Position of triangle vertex A.
  /// \param[in]  b Position of triangle vertex B.
  /// \param[in]  c Position of triangle vertex C.
  /// \param[out] n Non-normalized triangle normal vector. Length equals twice the triangle area.
  static void NormalDirection(const double a[3], const double b[3], const double c[3], double n[3]);

  /// Compute normal of triangle
  ///
  /// \param[in]  a Position of triangle vertex A.
  /// \param[in]  b Position of triangle vertex B.
  /// \param[in]  c Position of triangle vertex C.
  /// \param[out] n Triangle normal vector.
  static void Normal(const double a[3], const double b[3], const double c[3], double n[3]);

  /// Compute twice the area of triangle
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Twice the area of the triangle.
  static double DoubleArea(const double a[3], const double b[3], const double c[3]);

  /// Compute area of triangle
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Area of triangle.
  static double Area(const double a[3], const double b[3], const double c[3]);

  /// Compute area of triangle in 2D
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Signed area of triangle.
  static double SignedArea2D(const double a[2], const double b[2], const double c[2]);

  /// Compute twice the area of triangle in 2D
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Twice the signed area of the triangle.
  static double DoubleSignedArea2D(const double a[2], const double b[2], const double c[2]);

  /// Compute area of triangle in 2D
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Area of triangle.
  static double Area2D(const double a[2], const double b[2], const double c[2]);

  /// Compute twice the area of triangle in 2D
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Twice the area of the triangle.
  static double DoubleArea2D(const double a[2], const double b[2], const double c[2]);

  /// Compute cotangent of angle ABC
  ///
  /// \param[in] a Position of triangle vertex A.
  /// \param[in] b Position of triangle vertex B.
  /// \param[in] c Position of triangle vertex C.
  ///
  /// \returns Cotangent of angle ABC (equals cotangent of angle CBA).
  static double Cotangent(double a[3], double b[3], double c[3]);

  /// Tests whether two triangles intersect each other
  ///
  /// "A fast triangle to triangle intersection test for collision detection"
  /// Oren Tropp, Ayellet Tal, Ilan Shimshoni
  /// Computer Animation and Virtual Worlds 17(5) 2006, pp 527-535.
  static int TriangleTriangleIntersection(const double a1[3], const double b1[3], const double c1[3],
                                          const double a2[3], const double b2[3], const double c2[3]);

  /// Compute distance between closest points of two triangles
  static double DistanceBetweenTriangles(double a1[3], double b1[3], double c1[3], double n1[3],
                                         double a2[3], double b2[3], double c2[3], double n2[3],
                                         double *p1, double *p2);

  /// Compute distance between closest points of two triangles
  static double DistanceBetweenTriangles(const double a1[3], const double b1[3], const double c1[3], const double n1[3],
                                         const double a2[3], const double b2[3], const double c2[3], const double n2[3],
                                         double *p1, double *p2);

  /// Compute distance between closest points of two triangles
  static double DistanceBetweenTriangles(const double a1[3], const double b1[3], const double c1[3],
                                         const double a2[3], const double b2[3], const double c2[3],
                                         double *p1, double *p2);

};


////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// Copy 3D point coordinates stored in plain C array
#define MIRTK_RETURN_POINT(a, b) if (a) memcpy((a), (b), 3 * sizeof(double))

// =============================================================================
// Normals
// =============================================================================

// -----------------------------------------------------------------------------
inline void Triangle::NormalDirection(const double a[3], const double b[3], const double c[3], double n[3])
{
  double ab[3], ac[3];
  vtkMath::Subtract(b, a, ab);
  vtkMath::Subtract(c, a, ac);
  vtkMath::Cross(ab, ac, n);
}

// -----------------------------------------------------------------------------
inline void Triangle::Normal(const double a[3], const double b[3], const double c[3], double n[3])
{
  double l;
  NormalDirection(a, b, c, n);
  if ((l = vtkMath::Norm(n)) != 0.) {
    n[0] /= l, n[1] /= l, n[2] /= l;
  }
}

// =============================================================================
// Area
// =============================================================================

// -----------------------------------------------------------------------------
inline double Triangle::DoubleArea(const double a[3], const double b[3], const double c[3])
{
  double n[3];
  NormalDirection(a, b, c, n);
  return vtkMath::Norm(n);
}

// -----------------------------------------------------------------------------
inline double Triangle::Area(const double a[3], const double b[3], const double c[3])
{
  return .5 * DoubleArea(a, b, c);
}

// -----------------------------------------------------------------------------
inline double Triangle::DoubleSignedArea2D(const double a[2], const double b[2], const double c[2])
{
  return (a[0]*b[1] - a[1]*b[0]) + (b[0]*c[1] - b[1]*c[0]) + (c[0]*a[1] - c[1]*a[0]);
}

// -----------------------------------------------------------------------------
inline double Triangle::SignedArea2D(const double a[2], const double b[2], const double c[2])
{
  return .5 * DoubleSignedArea2D(a, b, c);
}

// -----------------------------------------------------------------------------
inline double Triangle::DoubleArea2D(const double a[2], const double b[2], const double c[2])
{
  return abs(DoubleSignedArea2D(a, b, c));
}

// -----------------------------------------------------------------------------
inline double Triangle::Area2D(const double a[2], const double b[2], const double c[2])
{
  return abs(SignedArea2D(a, b, c));
}

// =============================================================================
// Angles
// =============================================================================

// -----------------------------------------------------------------------------
// Meyer et al. (2002). Generalized Barycentric Coordinates on Irregular Polygons.
inline double Triangle::Cotangent(double a[3], double b[3], double c[3])
{
  double ba[3], bc[3], n[3];
  vtkMath::Subtract(a, b, ba);
  vtkMath::Subtract(c, b, bc);
  vtkMath::Cross(ba, bc, n);
  return vtkMath::Dot(ba, bc) / vtkMath::Norm(n);
}

// =============================================================================
// Triangle/triangle intersection test
// =============================================================================

// -----------------------------------------------------------------------------
/// Vector cross product
inline void CROSS(double *dest, const double *v1, const double *v2)
{
  dest[0] = v1[1] * v2[2] - v1[2] * v2[1];
  dest[1] = v1[2] * v2[0] - v1[0] * v2[2];
  dest[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

// -----------------------------------------------------------------------------
/// Vector dot product
inline double DOT(const double *v1, const double *v2)
{
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// -----------------------------------------------------------------------------
/// Subtract two vectors
inline void SUB(double *dest, const double *v1, const double *v2)
{
  dest[0] = v1[0] - v2[0];
  dest[1] = v1[1] - v2[1];
  dest[2] = v1[2] - v2[2];
}

// -----------------------------------------------------------------------------
// Scaled 2D vector addition
inline void sVpsV_2(double *Vr, double s1, const double *V1, double s2, const double *V2)
{
  Vr[0] = s1*V1[0] + s2*V2[0];
  Vr[1] = s1*V1[1] + s2*V2[1];
}

// -----------------------------------------------------------------------------
// Scaled vector subtraction
inline void VmsV(double *Vr, const double *V1, double s2, const double *V2)
{
  Vr[0] = V1[0] - s2 * V2[0];
  Vr[1] = V1[1] - s2 * V2[1];
  Vr[2] = V1[2] - s2 * V2[2];
}

// -----------------------------------------------------------------------------
// Sort numbers so that a<=b
inline void SORT(double &a, double &b)
{
  if( a > b) swap(a, b);
}

// -----------------------------------------------------------------------------
// Tomas Moller: this edge to edge test is based on Franlin Antonio's gem:
// "Faster Line Segment Intersection", in Graphics Gems III, pp. 199-202
#define MIRTK_EDGE_EDGE_TEST(V0,U0,U1)                                         \
  Bx=U0[i0]-U1[i0];                                                            \
  By=U0[i1]-U1[i1];                                                            \
  Cx=V0[i0]-U0[i0];                                                            \
  Cy=V0[i1]-U0[i1];                                                            \
  f=Ay*Bx-Ax*By;                                                               \
  d=By*Cx-Bx*Cy;                                                               \
  if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))                           \
  {                                                                            \
    e=Ax*Cy-Ay*Cx;                                                             \
    if(f>0)                                                                    \
    {                                                                          \
      if(e>=0 && e<=f) return 1;                                               \
    }                                                                          \
    else                                                                       \
    {                                                                          \
      if(e<=0 && e>=f) return 1;                                               \
    }                                                                          \
  }                                

// -----------------------------------------------------------------------------
#define MIRTK_EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2)                           \
{                                                                              \
  double Ax,Ay,Bx,By,Cx,Cy,e,d,f;                                              \
  Ax=V1[i0]-V0[i0];                                                            \
  Ay=V1[i1]-V0[i1];                                                            \
  /* test edge U0,U1 against V0,V1 */                                          \
  MIRTK_EDGE_EDGE_TEST(V0,U0,U1);                                              \
  /* test edge U1,U2 against V0,V1 */                                          \
  MIRTK_EDGE_EDGE_TEST(V0,U1,U2);                                              \
  /* test edge U2,U1 against V0,V1 */                                          \
  MIRTK_EDGE_EDGE_TEST(V0,U2,U0);                                              \
}

// -----------------------------------------------------------------------------
#define MIRTK_POINT_IN_TRI(V0,U0,U1,U2)                                        \
{                                                                              \
  double a,b,c,d0,d1,d2;                                                       \
  /* is T1 completly inside T2? */                                             \
  /* check if V0 is inside tri(U0,U1,U2) */                                    \
  a=U1[i1]-U0[i1];                                                             \
  b=-(U1[i0]-U0[i0]);                                                          \
  c=-a*U0[i0]-b*U0[i1];                                                        \
  d0=a*V0[i0]+b*V0[i1]+c;                                                      \
                                                                               \
  a=U2[i1]-U1[i1];                                                             \
  b=-(U2[i0]-U1[i0]);                                                          \
  c=-a*U1[i0]-b*U1[i1];                                                        \
  d1=a*V0[i0]+b*V0[i1]+c;                                                      \
                                                                               \
  a=U0[i1]-U2[i1];                                                             \
  b=-(U0[i0]-U2[i0]);                                                          \
  c=-a*U2[i0]-b*U2[i1];                                                        \
  d2=a*V0[i0]+b*V0[i1]+c;                                                      \
  if(d0*d1>0.0)                                                                \
  {                                                                            \
    if(d0*d2>0.0) return 1;                                                    \
  }                                                                            \
}

// -----------------------------------------------------------------------------
/// Procedure testing for intersection between coplanar triangles
///
/// Tomas Moller, "A Fast Triangle-Triangle Intersection Test",
/// Journal of Graphics Tools, 2(2), 1997
inline int coplanar_tri_tri(const double N[3],
                            const double V0[3], const double V1[3], const double V2[3],
                            const double U0[3], const double U1[3], const double U2[3])
{
  double A[3];
  short i0,i1;

  // first project onto an axis-aligned plane, that maximizes the area
  // of the triangles, compute indices: i0,i1.
  A[0]=abs(N[0]);
  A[1]=abs(N[1]);
  A[2]=abs(N[2]);
  if(A[0]>A[1]) {
    if(A[0]>A[2]) {
      i0=1; // A[0] is greatest
      i1=2;
    } else {
      i0=0; // A[2] is greatest
      i1=1;
    }
  } else {  // A[0]<=A[1]
    if(A[2]>A[1]) {
      i0=0; // A[2] is greatest
      i1=1;
    } else {
      i0=0; // A[1] is greatest
      i1=2;
    }
  }

  // test all edges of triangle 1 against the edges of triangle 2
  MIRTK_EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2);
  MIRTK_EDGE_AGAINST_TRI_EDGES(V1,V2,U0,U1,U2);
  MIRTK_EDGE_AGAINST_TRI_EDGES(V2,V0,U0,U1,U2);

  // finally, test if tri1 is totally contained in tri2 or vice versa
  MIRTK_POINT_IN_TRI(V0,U0,U1,U2);
  MIRTK_POINT_IN_TRI(U0,V0,V1,V2);

  return 0;
}

// -----------------------------------------------------------------------------
#define MIRTK_NEWCOMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,A,B,C,X0,X1) \
{                                                                              \
  if (D0D1>.0) {                                                               \
    /* here we know that D0D2<=0.0 */                                          \
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
    A=VV2; B=(VV0-VV2)*D2; C=(VV1-VV2)*D2; X0=D2-D0; X1=D2-D1;                 \
  } else if (D0D2>.0) {                                                        \
    /* here we know that d0d1<=0.0 */                                          \
    A=VV1; B=(VV0-VV1)*D1; C=(VV2-VV1)*D1; X0=D1-D0; X1=D1-D2;                 \
  } else if (D1*D2>.0 || D0!=.0) {                                             \
    /* here we know that d0d1<=0.0 or that D0!=0.0 */                          \
    A=VV0; B=(VV1-VV0)*D0; C=(VV2-VV0)*D0; X0=D0-D1; X1=D0-D2;                 \
  } else if (D1!=.0) {                                                         \
    A=VV1; B=(VV0-VV1)*D1; C=(VV2-VV1)*D1; X0=D1-D0; X1=D1-D2;                 \
  } else if (D2!=.0) {                                                         \
    A=VV2; B=(VV0-VV2)*D2; C=(VV1-VV2)*D2; X0=D2-D0; X1=D2-D1;                 \
  } else {                                                                     \
    /* triangles are coplanar */                                               \
    return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);                             \
  }                                                                            \
}

// -----------------------------------------------------------------------------
// Triangle/triangle intersection test routine by Tomas Moller, 1997.
// See article "A Fast Triangle-Triangle Intersection Test",
// Journal of Graphics Tools, 2(2), 1997
#define MIRTK_NoDivTriTriIsect_USE_EPSILON_TEST 0
inline int NoDivTriTriIsect(const double V0[3], const double V1[3], const double V2[3],
                            const double U0[3], const double U1[3], const double U2[3])
{
#if MIRTK_NoDivTriTriIsect_USE_EPSILON_TEST
  const double EPSILON = 0.000001;
#endif

  double E1[3],E2[3];
  double N1[3],N2[3],d1,d2;
  double du0,du1,du2,dv0,dv1,dv2;
  double D[3];
  double isect1[2], isect2[2];
  double du0du1,du0du2,dv0dv1,dv0dv2;
  short  index;
  double vp0,vp1,vp2;
  double up0,up1,up2;
  double bb,cc,max;

  /* compute plane equation of triangle(V0,V1,V2) */
  SUB(E1,V1,V0);
  SUB(E2,V2,V0);
  CROSS(N1,E1,E2);
  d1=-DOT(N1,V0);
  /* plane equation 1: N1.X+d1=0 */

  /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
  du0=DOT(N1,U0)+d1;
  du1=DOT(N1,U1)+d1;
  du2=DOT(N1,U2)+d1;

  /* coplanarity robustness check */
#if MIRTK_NoDivTriTriIsect_USE_EPSILON_TEST
  if(FABS(du0)<EPSILON) du0=.0;
  if(FABS(du1)<EPSILON) du1=.0;
  if(FABS(du2)<EPSILON) du2=.0;
#endif
  du0du1=du0*du1;
  du0du2=du0*du2;

  if(du0du1>.0 && du0du2>.0) /* same sign on all of them + not equal 0 ? */
    return 0;                /* no intersection occurs */

  /* compute plane of triangle (U0,U1,U2) */
  SUB(E1,U1,U0);
  SUB(E2,U2,U0);
  CROSS(N2,E1,E2);
  d2=-DOT(N2,U0);
  /* plane equation 2: N2.X+d2=0 */

  /* put V0,V1,V2 into plane equation 2 */
  dv0=DOT(N2,V0)+d2;
  dv1=DOT(N2,V1)+d2;
  dv2=DOT(N2,V2)+d2;

#if MIRTK_NoDivTriTriIsect_USE_EPSILON_TEST
  if(abs(dv0)<EPSILON) dv0=.0;
  if(abs(dv1)<EPSILON) dv1=.0;
  if(abs(dv2)<EPSILON) dv2=.0;
#endif

  dv0dv1=dv0*dv1;
  dv0dv2=dv0*dv2;

  if(dv0dv1>.0 && dv0dv2>.0) /* same sign on all of them + not equal 0 ? */
    return 0;                  /* no intersection occurs */

  /* compute direction of intersection line */
  CROSS(D,N1,N2);

  /* compute and index to the largest component of D */
  max=abs(D[0]);
  index=0;
  bb=abs(D[1]);
  cc=abs(D[2]);
  if(bb>max) max=bb,index=1;
  if(cc>max) max=cc,index=2;

  /* this is the simplified projection onto L*/
  vp0=V0[index];
  vp1=V1[index];
  vp2=V2[index];

  up0=U0[index];
  up1=U1[index];
  up2=U2[index];

  /* compute interval for triangle 1 */
  double a,b,c,x0,x1;
  MIRTK_NEWCOMPUTE_INTERVALS(vp0,vp1,vp2,dv0,dv1,dv2,dv0dv1,dv0dv2,a,b,c,x0,x1);

  /* compute interval for triangle 2 */
  double d,e,f,y0,y1;
  MIRTK_NEWCOMPUTE_INTERVALS(up0,up1,up2,du0,du1,du2,du0du1,du0du2,d,e,f,y0,y1);

  double xx,yy,xxyy,tmp;
  xx=x0*x1;
  yy=y0*y1;
  xxyy=xx*yy;

  tmp=a*xxyy;
  isect1[0]=tmp+b*x1*yy;
  isect1[1]=tmp+c*x0*yy;

  tmp=d*xxyy;
  isect2[0]=tmp+e*xx*y1;
  isect2[1]=tmp+f*xx*y0;

  SORT(isect1[0],isect1[1]);
  SORT(isect2[0],isect2[1]);


  if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return 0;
  return 1;
}
#undef MIRTK_NoDivTriTriIsect_USE_EPSILON_TEST

// -----------------------------------------------------------------------------
/// "A fast triangle to triangle intersection test for collision detection"
///
/// Oren Tropp, Ayellet Tal, Ilan Shimshoni
/// Computer Animation and Virtual Worlds 17(5) 2006, pp 527-535.
///
/// C1,C2,C3 are the vertices of triangle A. P1,P2 are the two edges originating from C1.
/// D1,D2,D3 are the vertices of triangle B. Q1,Q2 are the two edges originating from D1.
///
/// @return Zero for disjoint triangles and non-zero for intersection.
inline int tri_tri_intersect3D(const double *C1, const double *C2, const double *C3,
                               const double *P1, const double *P2,
	                             const double *D1, const double *D2, const double *D3,
                               const double *Q1, const double *Q2)
{
	double t[3],p1[3], p2[3],r[3],r4[3];
	double beta1, beta2, beta3;
	double gama1, gama2, gama3;
	double det1, det2, det3;
	double dp0, dp1, dp2;
	double dq1,dq2,dq3,dr, dr3;
	double alpha1, alpha2;
	bool   alpha1_legal, alpha2_legal;
	double SF;
	bool   beta1_legal, beta2_legal;

	SUB(r, D1, C1);

	// determinant computation	
	dp0 = P1[1]*P2[2]-P2[1]*P1[2];
	dp1 = P1[0]*P2[2]-P2[0]*P1[2];
	dp2 = P1[0]*P2[1]-P2[0]*P1[1];
	dq1 = Q1[0]*dp0 - Q1[1]*dp1 + Q1[2]*dp2;
	dq2 = Q2[0]*dp0 - Q2[1]*dp1 + Q2[2]*dp2;
	dr  = -r[0]*dp0  + r[1]*dp1  - r[2]*dp2;

	beta1 = dr*dq2;  // beta1, beta2 are scaled so that beta_i=beta_i*dq1*dq2
	beta2 = dr*dq1;
	beta1_legal = (beta2>=0) && (beta2 <=dq1*dq1) && (dq1 != 0);
	beta2_legal = (beta1>=0) && (beta1 <=dq2*dq2) && (dq2 != 0);
		
	dq3=dq2-dq1;
	dr3=+dr-dq1;   // actually this is -dr3

	if ((dq1 == .0) && (dq2 == .0))
	{
		if (dr != .0) return 0;  // triangles are on parallel planes
		// triangles are on the same plane, use the coplanar test of Moller
    double N1[3];
    CROSS(N1, P1, P2);
    return coplanar_tri_tri(N1, C1,C2,C3, D1,D2,D3);
	}

	if (!beta2_legal && !beta1_legal) return 0; // fast reject-all vertices are on
													                    // the same side of the triangle plane

	if (beta2_legal && beta1_legal) {   //beta1, beta2
		SF = dq1*dq2;
		sVpsV_2(t, beta2, Q2, (-beta1), Q1);
	} else if (beta1_legal && !beta2_legal) {   //beta1, beta3
		SF    = dq1*dq3;
		beta1 = beta1-beta2;   // all betas are multiplied by a positive SF
		beta3 = dr3*dq1;
		sVpsV_2(t, (SF-beta3-beta1), Q1, beta3, Q2);
	} else {   // beta2, beta3
		SF    = dq2*dq3;
		beta2 = beta1-beta2;   // all betas are multiplied by a positive SF
		beta3 = dr3*dq2;
		sVpsV_2(t, (SF-beta3), Q1, (beta3-beta2), Q2);
		Q1    = Q2;
		beta1 = beta2;
	}
	sVpsV_2(r4, SF, r, beta1, Q1);

  // seg_collide3(t, r4); // calculates the 2D intersection
  p1[0] = SF*P1[0];
  p1[1] = SF*P1[1];
  p2[0] = SF*P2[0];
  p2[1] = SF*P2[1];

  det1         = p1[0]*t[1]-t[0]*p1[1];
  gama1        = (p1[0]*r4[1]-r4[0]*p1[1])*det1;
  alpha1       = (r4[0]*t[1] - t[0]*r4[1])*det1;
  alpha1_legal = (alpha1>=0) && (alpha1<=(det1*det1)  && (det1!=0));
  det2         = p2[0]*t[1] - t[0]*p2[1];
  alpha2       = (r4[0]*t[1] - t[0]*r4[1]) *det2;
  gama2        = (p2[0]*r4[1] - r4[0]*p2[1]) * det2;
  alpha2_legal = (alpha2>=0) && (alpha2<=(det2*det2) && (det2 !=0));
  det3         = det2-det1;
  gama3        = ((p2[0]-p1[0])*(r4[1]-p1[1]) - (r4[0]-p1[0])*(p2[1]-p1[1]))*det3;

  if (alpha1_legal) {
    if (alpha2_legal) {
      if ( ((gama1<=0) && (gama1>=-(det1*det1))) ||
           ((gama2<=0) && (gama2>=-(det2*det2))) || (gama1*gama2<0)) return 12;
    } else {
      if ( ((gama1<=0) && (gama1>=-(det1*det1))) ||
           ((gama3<=0) && (gama3>=-(det3*det3))) || (gama1*gama3<0)) return 13;
    }
  } else {
    if (alpha2_legal) {
      if ( ((gama2<=0) && (gama2>=-(det2*det2))) ||
           ((gama3<=0) && (gama3>=-(det3*det3))) || (gama2*gama3<0)) return 23;
    }
  }
	return 0;
}

// -----------------------------------------------------------------------------
inline int Triangle
::TriangleTriangleIntersection(const double c1[3], const double c2[3], const double c3[3],
                               const double d1[3], const double d2[3], const double d3[3])
{
  const double p1[3] = { c2[0] - c1[0], c2[1] - c1[1], c2[2] - c1[2] };
  const double p2[3] = { c3[0] - c1[0], c3[1] - c1[1], c3[2] - c1[2] };
  const double q1[3] = { d2[0] - d1[0], d2[1] - d1[1], d2[2] - d1[2] };
  const double q2[3] = { d3[0] - d1[0], d3[1] - d1[1], d3[2] - d1[2] };
  return tri_tri_intersect3D(c1, c2, c3, p1, p2, d1, d2, d3, q1, q2);
}

// =============================================================================
// Distance between two triangles
// =============================================================================

// -----------------------------------------------------------------------------
/// Compute distance between two triangles
inline double Triangle
::DistanceBetweenTriangles(double a1[3], double b1[3], double c1[3], double n1[3],
                           double a2[3], double b2[3], double c2[3], double n2[3],
                           double *p1, double *p2)
{
  const double TOL = 1e-6; // Tolerance for point coordinate equality check
  double q1[3], q2[3], t1, t2, d, min = numeric_limits<double>::infinity();

  // Point in second triangle closest to any of the corners of the first triangle
  t2 = vtkPlane::Evaluate(n2, a2, a1);
  d  = abs(t2);
  if (d < min) {
    VmsV(q2, a1, t2, n2);
    if (vtkTriangle::PointInTriangle(q2, a2, b2, c2, TOL)) {
      MIRTK_RETURN_POINT(p1, a1);
      MIRTK_RETURN_POINT(p2, q2);
      min = d;
    }
  }
  t2 = vtkPlane::Evaluate(n2, a2, b1);
  d  = abs(t2);
  if (d < min) {
    VmsV(q2, b1, t2, n2);
    if (vtkTriangle::PointInTriangle(q2, a2, b2, c2, TOL)) {
      MIRTK_RETURN_POINT(p1, b1);
      MIRTK_RETURN_POINT(p2, q2);
      min = d;
    }
  }
  t2 = vtkPlane::Evaluate(n2, a2, c1);
  d  = abs(t2);
  if (d < min) {
    VmsV(q2, c1, t2, n2);
    if (vtkTriangle::PointInTriangle(q2, a2, b2, c2, TOL)) {
      MIRTK_RETURN_POINT(p1, c1);
      MIRTK_RETURN_POINT(p2, q2);
      min = d;
    }
  }

  // Point in first triangle closest to any of the corners of the second triangle
  t1 = vtkPlane::Evaluate(n1, a1, a2);
  d  = abs(t1);
  if (d < min) {
    VmsV(q1, a2, t1, n1);
    if (vtkTriangle::PointInTriangle(q1, a1, b1, c1, TOL)) {
      MIRTK_RETURN_POINT(p1, q1);
      MIRTK_RETURN_POINT(p2, a2);
      min = d;
    }
  }
  t1 = vtkPlane::Evaluate(n1, a1, b2);
  d  = abs(t1);
  if (d < min) {
    VmsV(q1, b2, t1, n1);
    if (vtkTriangle::PointInTriangle(q1, a1, b1, c1, TOL)) {
      MIRTK_RETURN_POINT(p1, q1);
      MIRTK_RETURN_POINT(p2, b2);
      min = d;
    }
  }
  t1 = vtkPlane::Evaluate(n1, a1, c2);
  d  = abs(t1);
  if (d < min) {
    VmsV(q1, c2, t1, n1);
    if (vtkTriangle::PointInTriangle(q1, a1, b1, c1, TOL)) {
      MIRTK_RETURN_POINT(p1, q1);
      MIRTK_RETURN_POINT(p2, c2);
      min = d;
    }
  }

  // Distance between each pair of edges
  d = vtkLine::DistanceBetweenLineSegments(a1, b1, a2, b2, q1, q2, t1, t2);
  if (d < min) {
    MIRTK_RETURN_POINT(p1, q1);
    MIRTK_RETURN_POINT(p2, q2);
    min = d;
  }
  d = vtkLine::DistanceBetweenLineSegments(a1, b1, b2, c2, q1, q2, t1, t2);
  if (d < min) {
    MIRTK_RETURN_POINT(p1, q1);
    MIRTK_RETURN_POINT(p2, q2);
    min = d;
  }
  d = vtkLine::DistanceBetweenLineSegments(a1, b1, c2, a2, q1, q2, t1, t2);
  if (d < min) {
    MIRTK_RETURN_POINT(p1, q1);
    MIRTK_RETURN_POINT(p2, q2);
    min = d;
  }
  d = vtkLine::DistanceBetweenLineSegments(b1, c1, a2, b2, q1, q2, t1, t2);
  if (d < min) {
    MIRTK_RETURN_POINT(p1, q1);
    MIRTK_RETURN_POINT(p2, q2);
    min = d;
  }
  d = vtkLine::DistanceBetweenLineSegments(b1, c1, b2, c2, q1, q2, t1, t2);
  if (d < min) {
    MIRTK_RETURN_POINT(p1, q1);
    MIRTK_RETURN_POINT(p2, q2);
    min = d;
  }
  d = vtkLine::DistanceBetweenLineSegments(b1, c1, c2, a2, q1, q2, t1, t2);
  if (d < min) {
    MIRTK_RETURN_POINT(p1, q1);
    MIRTK_RETURN_POINT(p2, q2);
    min = d;
  }
  d = vtkLine::DistanceBetweenLineSegments(c1, a1, a2, b2, q1, q2, t1, t2);
  if (d < min) {
    MIRTK_RETURN_POINT(p1, q1);
    MIRTK_RETURN_POINT(p2, q2);
    min = d;
  }
  d = vtkLine::DistanceBetweenLineSegments(c1, a1, b2, c2, q1, q2, t1, t2);
  if (d < min) {
    MIRTK_RETURN_POINT(p1, q1);
    MIRTK_RETURN_POINT(p2, q2);
    min = d;
  }
  d = vtkLine::DistanceBetweenLineSegments(c1, a1, c2, a2, q1, q2, t1, t2);
  if (d < min) {
    MIRTK_RETURN_POINT(p1, q1);
    MIRTK_RETURN_POINT(p2, q2);
    min = d;
  }

  return min;
}

// -----------------------------------------------------------------------------
/// Compute distance between two triangles
inline double Triangle
::DistanceBetweenTriangles(const double a1[3], const double b1[3],
                           const double c1[3], const double n1[3],
                           const double a2[3], const double b2[3],
                           const double c2[3], const double n2[3],
                           double      *p1,    double      *p2)
{
  return DistanceBetweenTriangles(const_cast<double *>(a1),
                                  const_cast<double *>(b1),
                                  const_cast<double *>(c1),
                                  const_cast<double *>(n1),
                                  const_cast<double *>(a2),
                                  const_cast<double *>(b2),
                                  const_cast<double *>(c2),
                                  const_cast<double *>(n2), p1, p2);
}

// -----------------------------------------------------------------------------
/// Compute distance between two triangles
inline double Triangle
::DistanceBetweenTriangles(const double a1[3], const double b1[3], const double c1[3],
                           const double a2[3], const double b2[3], const double c2[3],
                           double *p1, double *p2)
{
  double n1[3], n2[3];
  Normal(a1, b1, c1, n1);
  Normal(a2, b2, c2, n2);
  return DistanceBetweenTriangles(a1, b1, c1, n1, a2, b2, c2, n2, p1, p2);
}


// -----------------------------------------------------------------------------
// Undefine auxiliary macros again
#undef MIRTK_RETURN_POINT
#undef MIRTK_EDGE_EDGE_TEST
#undef MIRTK_EDGE_AGAINST_TRI_EDGES
#undef MIRTK_POINT_IN_TRI
#undef MIRTK_NEWCOMPUTE_INTERVALS


} // namespace mirtk

#endif // MIRTK_Triangle_H
