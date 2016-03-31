/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
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

#include "mirtk/EuclideanDistanceTransform.h"

#include "mirtk/Math.h"
#include "mirtk/Stream.h"

#define EDT_MAX_IMAGE_DIMENSION 26754
#define EDT_MAX_DISTANCE_SQUARED 2147329548
#define EDT_MAX_DISTANCE_SQUARED_ANISOTROPIC 2147329548


namespace mirtk {


// -----------------------------------------------------------------------------
template <class VoxelType>
EuclideanDistanceTransform<VoxelType>
::EuclideanDistanceTransform(Mode distanceTransformMode)
{
  _distanceTransformMode = distanceTransformMode;
}

// -----------------------------------------------------------------------------
// This procedure computes the squared EDT of a 2D binary image with isotropic
// voxels of unit dimension. The EDT can obviously be obtained simply from the
// output by taking the square root of each element. The output can be scaled
// to account for non-unit dimension. But neither of these functions are
// provided as options in this procedure.
//
// The EDT is returned in the array edt. Memory for the array edt must be
// allocated by the caller. The binary image can be provided in two different
// ways. It can be provided in the array img. In this case it is not changed by
// the procedure. Alternatively, the binary image can be provided in the array
// edt. In this case, the array img must be the NULL pointer, and the binary
// image will be overwritten by the EDT. The binary image doesn't have to
// consist of 0's and 1's. A voxel value of 0 denotes a background voxel and
// any other value denotes a foreground or feature voxel. The binary image and
// EDT dimensions are nX x nY voxels.
//
// The procedure uses the algorithm described in the paper:
// CR Maurer Jr, R Qi, V Raghavan. A linear time algorithm for computing exact
// Euclidean distance transforms of binary images in arbitrary dimensions.
// IEEE Transactions on Pattern Analysis and Machine Intelligence. In review.
// A preliminary version of this paper was published in the conference
// proceedings:
// CR Maurer Jr, V Raghavan, R Qi. A linear time algorithm for computing the
// Euclidean distance transform in arbitrary dimensions. In: MF Insana, RM
// Leahy, eds. Information Processing in Medical Imaging (IPMI) 2001. Berlin:
// Springer-Verlag, 2001, pp. 358-364. (Davis, CA, June 18-22, 2001).
template <class VoxelType>
void EuclideanDistanceTransform<VoxelType>::edtComputeEDT_2D(char *img, long *edt, long nX, long nY)
{
  const long max_dist2 = long(EDT_MAX_DISTANCE_SQUARED);

  char *c;
  long i, j, nXY, d, *p, *q, *f;

  /* nXY is number of voxels in 2D image */
  nXY = nX * nY;

  /* if binary image is provided in the array img, copy it to the arry edt */
  /* this is effectively equivalent to computing D_0 */
  if (img != NULL) {
    c = img;
    p = edt;
    for (i = 0; i < nXY; i++, c++, p++) {
      *p = *c;
    }
  }

  /* compute D_1 as simple forward-and-reverse distance propagation */
  /* (instead of calling edtVornoiEDT) */
  /* D_1 is distance to closest feature voxel in row (x direction) */
  /* it is possible to use a simple distance propagation for D_1  because */
  /* L_1 and L_2 norms are equivalent for 1D case */
  for (j = 0; j < nY; j++) {
    /* forward pass */
    p = edt + j * nX;
    d = max_dist2;
    for (i = 0; i < nX; i++, p++) {
      /* set d = 0 when we encounter a feature voxel */
      if (*p) {
        *p = d = 0;
      }
      /* increment distance ... */
      else if (d != max_dist2) {
        *p = ++d;
      }
      /* ... unless we haven't encountered a feature voxel yet */
      else {
        *p = max_dist2;
      }
    }
    /* reverse pass */
    if (*(--p) != max_dist2) {
      d = max_dist2;
      for (i = nX - 1; i >= 0; i--, p--) {
        /* set d = 0 when we encounter a feature voxel */
        if (*p == 0) {
          d = 0;
        }
        /* increment distance after encountering a feature voxel */
        else if (d != max_dist2) {
          /* compare forward and reverse distances */
          if (++d < *p) {
            *p = d;
          }
        }
        /* square distance */
        /* (we use squared distance in rest of algorithm) */
        *p *= *p;
      }
    }
  }

  /* compute D_2 = squared EDT */
  /* solve 1D problem for each column (y direction) */
  f = (long *)malloc(nY * sizeof(long));
  if (f == NULL) {
    fprintf(stderr, "Error in edtComputeEDT_2D()\n");
    fprintf(stderr, "Cannot malloc f\n");
    exit(EXIT_FAILURE);
  }
  for (i = 0; i < nX; i++) {
    /* fill array f with D_1 distances in column */
    /* this is essentially line 4 in Procedure VoronoiEDT() in tPAMI paper */
    p = edt + i;
    q = f;
    for (j = 0; j < nY; j++, p += nX, q++) {
      *q = *p;
    }
    /* call edtVornoiEDT */
    if (edtVornoiEDT(f, nY)) {
      p = edt + i;
      q = f;
      for (j = 0; j < nY; j++, p += nX, q++) {
        *p = *q;
      }
    }
  }
  free(f);
} /* edtComputeEDT_2D */

// -----------------------------------------------------------------------------
// This procedure computes the squared EDT of a 3D binary image with isotropic
// voxels of unit dimension. See notes for edtComputeEDT_2D.
template <class VoxelType>
void EuclideanDistanceTransform<VoxelType>::edtComputeEDT_3D(char *img, long *edt, long nX, long nY, long nZ)
{
  char *c;
  long i, k, nXY, nXYZ, *p, *q, *f;

  /* nXY is number of voxels in each plane (xy) */
  /* nXYZ is number of voxels in 3D image */
  nXY = nX * nY;
  nXYZ = nX * nY * nZ;

  /* if binary image is provided in the array img, copy it to the arry edt */
  /* this is effectively equivalent to computing D_0 */
  if (img != NULL) {
    c = img;
    p = edt;
    for (i = 0; i < nXYZ; i++, c++, p++) {
      *p = *c;
    }
  }

  /* compute D_2 */
  /* call edtComputeEDT_2D for each plane */
  p = edt;
  for (k = 0; k < nZ; k++, p += nXY) {
    edtComputeEDT_2D(NULL, p, nX, nY);
  }

  /* compute D_3 */
  /* solve 1D problem for each column (z direction) */
  f = (long *)malloc(nZ * sizeof(long));
  if (f == NULL) {
    fprintf(stderr, "Error in edtComputeEDT_3D()\n");
    fprintf(stderr, "Cannot malloc f\n");
    exit(EXIT_FAILURE);
  }
  for (i = 0; i < nXY; i++) {
    /* fill array f with D_2 distances in column */
    /* this is essentially line 4 in Procedure VoronoiEDT() in tPAMI paper */
    p = edt + i;
    q = f;
    for (k = 0; k < nZ; k++, p += nXY, q++) {
      *q = *p;
    }
    /* call edtVornoiEDT */
    if (edtVornoiEDT(f, nZ)) {
      p = edt + i;
      q = f;
      for (k = 0; k < nZ; k++, p += nXY, q++) {
        *p = *q;
      }
    }
  }
  free(f);
} /* edtComputeEDT_3D */

// -----------------------------------------------------------------------------
// This is Procedure edtVornoiEDT() in tPAMI paper.
template <class VoxelType>
int EuclideanDistanceTransform<VoxelType>::edtVornoiEDT(long *f, long n)
{
  const VoxelType max_dist2 = VoxelType(EDT_MAX_DISTANCE_SQUARED);

  long i, l, a, b, c, v, n_S, lhs, rhs;
  static int firstCall = 1;
  static long *g, *h, size_gh = 0;

  /* this procedure is called often */
  /* we keep static arrays to avoid frequent calls to malloc() */
  /* the downside is that this memory is never freed, */
  /* but it is not so much memory: max(nX,nY,nZ)*sizeof(long) */
  /* need to check if arrays are sufficiently large */
  if (n > size_gh && !firstCall) {
    free(g);
    free(h);
    firstCall = 1;
  }
  /* malloc arrays if this is first call to procedure, or if arrays */
  /* are too small and need to be reallocated */
  if (firstCall) {
    g = (long *)malloc(n * sizeof(long));
    if (g == NULL) {
      fprintf(stderr, "Error in edtVornoiEDT()\n");
      fprintf(stderr, "Cannot malloc g\n");
      exit(EXIT_FAILURE);
    }
    h = (long *)malloc(n * sizeof(long));
    if (h == NULL) {
      fprintf(stderr, "Error in edtVornoiEDT()\n");
      fprintf(stderr, "Cannot malloc h\n");
      exit(EXIT_FAILURE);
    }
    size_gh = n;
    firstCall = 0;
  }

  /* construct partial Vornoi diagram */
  /* this loop is lines 1-14 in Procedure edtVornoiEDT() in tPAMI paper */
  /* note we use 0 indexing in this program whereas paper uses 1 indexing */
  for (i = 0, l = -1; i < n; i++) {
    /* line 4 */
    if (f[i] != max_dist2) {
      /* line 5 */
      if (l < 1) {
        /* line 6 */
        g[++l] = f[i];
        h[l] = i;
      }
      /* line 7 */
      else {
        /* line 8 */
        while (l >= 1) {
          /* compute removeEDT() in line 8 */
          v = h[l];
          a = v - h[l-1];
          b = i - v;
          c = a + b;
          /* compute Eq. 2 */
          if ((c*g[l] - b*g[l-1] - a*f[i] - a*b*c) > 0) {
            /* line 9 */
            l--;
          } else {
            break;
          }
        }
        /* line 11 */
        g[++l] = f[i];
        h[l] = i;
      }
    }
  }
  /* query partial Vornoi diagram */
  /* this is lines 15-25 in Procedure edtVornoiEDT() in tPAMI paper */
  /* lines 15-17 */
  if ((n_S = l + 1) == 0) {
    return (0);
  }
  /* lines 18-19 */
  for (i = 0, l = 0; i < n; i++) {
    /* line 20 */
    /* we reduce number of arithmetic operations by taking advantage of */
    /* similarities in successive computations instead of treating them as */
    /* independent ones */
    a = h[l] - i;
    lhs = g[l] + a * a;
    while (l < n_S - 1) {
      a = h[l+1] - i;
      rhs = g[l+1] + a * a;
      if (lhs > rhs) {
        /* line 21 */
        l++;
        lhs = rhs;
      } else {
        break;
      }
    }
    /* line 23 */
    /* we put distance into the 1D array that was passed; */
    /* must copy into EDT in calling procedure */
    f[i] = lhs;
  }
  /* line 25 */
  /* return 1 if we queried diagram, 0 if we returned because n_S = 0 */
  return (1);
} /* edtVornoiEDT */

// -----------------------------------------------------------------------------
// This is Procedure edtVornoiEDT() in tPAMI paper.
template <class VoxelType>
int EuclideanDistanceTransform<VoxelType>::edtVornoiEDT_anisotropic(VoxelType *f, long n, double w)
{
  const VoxelType max_dist2 = VoxelType(EDT_MAX_DISTANCE_SQUARED_ANISOTROPIC);

  long i, l, n_S;
  float a, b, c, v, lhs, rhs;
  static int firstCall = 1;
  static long size_gh = 0;
  static float *g, *h;

  /* this procedure is called often */
  /* we keep static arrays to avoid frequent calls to malloc() */
  /* the downside is that this memory is never freed, */
  /* but it is not so much memory: max(nX,nY,nZ)*sizeof(long) */
  /* need to check if arrays are sufficiently large */
  if (n > size_gh && !firstCall) {
    free(g);
    free(h);
    firstCall = 1;
  }
  /* malloc arrays if this is first call to procedure, or if arrays */
  /* are too small and need to be reallocated */
  if (firstCall) {
    g = (float *)malloc(n * sizeof(float));
    if (g == NULL) {
      fprintf(stderr, "Error in edtVornoiEDT()\n");
      fprintf(stderr, "Cannot malloc g\n");
      exit(EXIT_FAILURE);
    }
    h = (float *)malloc(n * sizeof(float));
    if (h == NULL) {
      fprintf(stderr, "Error in edtVornoiEDT()\n");
      fprintf(stderr, "Cannot malloc h\n");
      exit(EXIT_FAILURE);
    }
    size_gh = n;
    firstCall = 0;
  }

  /* construct partial Vornoi diagram */
  /* this loop is lines 1-14 in Procedure edtVornoiEDT() in tPAMI paper */
  /* note we use 0 indexing in this program whereas paper uses 1 indexing */
  for (i = 0, l = -1; i < n; i++) {
    /* line 4 */
    if (f[i] != max_dist2) {
      /* line 5 */
      if (l < 1) {
        /* line 6 */
        g[++l] = static_cast<float>(f[i]);
        h[l] = static_cast<float>(w * i);
      }
      /* line 7 */
      else {
        /* line 8 */
        while (l >= 1) {
          /* compute removeEDT() in line 8 */
          v = h[l];
          a = v - h[l-1];
          b = static_cast<float>(w * i) - v;
          c = a + b;
          /* compute Eq. 2 */
          if ((c*g[l] - b*g[l-1] - a*static_cast<float>(f[i]) - a*b*c) > .0f) {
            /* line 9 */
            l--;
          } else {
            break;
          }
        }
        /* line 11 */
        g[++l] = static_cast<float>(f[i]);
        h[l] = static_cast<float>(w * i);
      }
    }
  }
  /* query partial Vornoi diagram */
  /* this is lines 15-25 in Procedure edtVornoiEDT() in tPAMI paper */
  /* lines 15-17 */
  if ((n_S = l + 1) == 0) {
    return (0);
  }
  /* lines 18-19 */
  for (i = 0, l = 0; i < n; i++) {
    /* line 20 */
    /* we reduce number of arithmetic operations by taking advantage of */
    /* similarities in successive computations instead of treating them as */
    /* independent ones */
    a = h[l] - static_cast<float>(w * i);
    lhs = g[l] + a * a;
    while (l < n_S - 1) {
      a = h[l+1] - static_cast<float>(w * i);
      rhs = g[l+1] + a * a;
      if (lhs > rhs) {
        /* line 21 */
        l++;
        lhs = rhs;
      } else {
        break;
      }
    }
    /* line 23 */
    /* we put distance into the 1D array that was passed; */
    /* must copy into EDT in calling procedure */
    f[i] = static_cast<VoxelType>(lhs);
  }
  /* line 25 */
  /* return 1 if we queried diagram, 0 if we returned because n_S = 0 */
  return (1);
} /* edtVornoiEDT_anisotropic */

// -----------------------------------------------------------------------------
// This procedure computes the squared EDT of a 2D binary image with anisotropic
// voxels. See notes for edtComputeEDT_2D. The difference relative to edtComputeEDT_2D
// is that the edt is a float array instead of a long array, and there are
// additional parameters for the image voxel dimensions wX and wY.
template <class VoxelType>
void EuclideanDistanceTransform<VoxelType>
::edtComputeEDT_2D_anisotropic(const VoxelType *img, VoxelType *edt, long nX, long nY, double wX, double wY)
{
  const VoxelType max_dist2 = VoxelType(EDT_MAX_DISTANCE_SQUARED_ANISOTROPIC);

  const VoxelType *c;
  VoxelType d, *p, *q, *f;

  /* nXY is number of voxels in 2D image */
  const long nXY = nX * nY;

  /* if binary image is provided in the array img, copy it to the arry edt */
  /* this is effectively equivalent to computing D_0 */
  if (img != nullptr) {
    c = img;
    p = edt;
    memcpy(edt, img, nXY * sizeof(VoxelType));
  }

  /* compute D_1 as simple forward-and-reverse distance propagation */
  /* (instead of calling edtVornoiEDT) */
  /* D_1 is distance to closest feature voxel in row (x direction) */
  /* it is possible to use a simple distance propagation for D_1  because */
  /* L_1 and L_2 norms are equivalent for 1D case */
  for (long j = 0; j < nY; j++) {
    /* forward pass */
    p = edt + j * nX;
    d = max_dist2;
    for (long i = 0; i < nX; i++, p++) {
      /* set d = 0 when we encounter a feature voxel */
      if (*p) {
        *p = d = 0;
      /* increment distance ... */
      } else if (d != max_dist2) {
        *p = ++d;
      /* ... unless we haven't encountered a feature voxel yet */
      } else {
        *p = max_dist2;
      }
    }
    /* reverse pass */
    if (*(--p) != max_dist2) {
      d = max_dist2;
      for (long i = nX - 1; i >= 0; i--, p--) {
        /* set d = 0 when we encounter a feature voxel */
        if (*p == 0) {
          d = 0;
        /* increment distance after encountering a feature voxel */
        } else if (d != max_dist2) {
          /* compare forward and reverse distances */
          if (++d < *p) {
            *p = d;
          }
        }
        /* square distance */
        /* (we use squared distance in rest of algorithm) */
        *p *= static_cast<VoxelType>(wX);
        *p *= *p;
      }
    }
  }

  /* compute D_2 = squared EDT */
  /* solve 1D problem for each column (y direction) */
  f = (VoxelType *)malloc(nY * sizeof(VoxelType));
  if (f == NULL) {
    fprintf(stderr, "Error in edtComputeEDT_2D()\n");
    fprintf(stderr, "Cannot malloc f\n");
    exit(EXIT_FAILURE);
  }
  for (long i = 0; i < nX; i++) {
    /* fill array f with D_1 distances in column */
    /* this is essentially line 4 in Procedure VoronoiEDT() in tPAMI paper */
    p = edt + i;
    q = f;
    for (long j = 0; j < nY; j++, p += nX, q++) {
      *q = *p;
    }
    /* call edtVornoiEDT */
    if (edtVornoiEDT_anisotropic(f, nY, wY)) {
      p = edt + i;
      q = f;
      for (long j = 0; j < nY; j++, p += nX, q++) {
        *p = *q;
      }
    }
  }
  free(f);
} /* edtComputeEDT_2D_anisotropic */

// -----------------------------------------------------------------------------
// This procedure computes the squared EDT of a 3D binary image with anisotropic
// voxels. See notes for edtComputeEDT_2D_anisotropic.
template <class VoxelType>
void EuclideanDistanceTransform<VoxelType>
::edtComputeEDT_3D_anisotropic(const VoxelType *img, VoxelType *edt,
                               long nX, long nY, long nZ, double wX, double wY, double wZ)
{
  const VoxelType *c;
  long i, k, nXY, nXYZ;
  VoxelType *p, *q, *f;

  /* nXY is number of voxels in each plane (xy) */
  /* nXYZ is number of voxels in 3D image */
  nXY = nX * nY;
  nXYZ = nX * nY * nZ;

  /* if binary image is provided in the array img, copy it to the arry edt */
  /* this is effectively equivalent to computing D_0 */
  if (img != NULL) {
    c = img;
    p = edt;
    for (i = 0; i < nXYZ; i++, c++, p++) {
      *p = *c;
    }
  }

  /* compute D_2 */
  /* call edtComputeEDT_2D for each plane */
  p = edt;
  for (k = 0; k < nZ; k++, p += nXY) {
    edtComputeEDT_2D_anisotropic(NULL, p, nX, nY, wX, wY);
  }

  /* compute D_3 */
  /* solve 1D problem for each column (z direction) */
  f = (VoxelType *)malloc(nZ * sizeof(VoxelType));
  if (f == NULL) {
    fprintf(stderr, "Error in edtComputeEDT_3D()\n");
    fprintf(stderr, "Cannot malloc f\n");
    exit(EXIT_FAILURE);
  }
  for (i = 0; i < nXY; i++) {
    /* fill array f with D_2 distances in column */
    /* this is essentially line 4 in Procedure VoronoiEDT() in tPAMI paper */
    p = edt + i;
    q = f;
    for (k = 0; k < nZ; k++, p += nXY, q++) {
      *q = *p;
    }
    /* call edtVornoiEDT */
    if (edtVornoiEDT_anisotropic(f, nZ, wZ)) {
      p = edt + i;
      q = f;
      for (k = 0; k < nZ; k++, p += nXY, q++) {
        *p = *q;
      }
    }
  }
  free(f);
} /* edtComputeEDT_3D_anisotropic */

// -----------------------------------------------------------------------------
template <class VoxelType>
void EuclideanDistanceTransform<VoxelType>::Run()
{
  // Do the initial set up
  this->Initialize();

  // Get image pointers (*after* Initialize()!)
  const ImageType * const input  = this->Input();
  ImageType       * const output = this->Output();

  // Get image size
  const int nx = input->X();
  const int ny = input->Y();
  const int nz = input->Z();
  const int nt = input->T();

  // Get voxel size
  const double dx = input->XSize();
  const double dy = input->YSize();
  const double dz = input->ZSize();

  for (int t = 0; t < nt; ++t) {
	  if (_distanceTransformMode == DT_3D) {
		  edtComputeEDT_3D_anisotropic(
        input ->Data(0, 0, 0, t),
			  output->Data(0, 0, 0, t),
			  nx, ny, nz, dx, dy, dz);
	  } else {
		  for (int z = 0; z < nz; ++z) {
			  edtComputeEDT_2D_anisotropic(
          input ->Data(0, 0, z, t),
				  output->Data(0, 0, z, t),
				  nx, ny, dx, dy);
		  }
	  }
  }

  // Do the final cleaning up
  this->Finalize();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void EuclideanDistanceTransform<VoxelType>::Radial()
{
  int nx, ny, nz, x, y, z;
  double min;

  // Calculate image dimensions
  nx = this->_Input->X();
  ny = this->_Input->Y();
  nz = this->_Input->Z();

  if (this->_distanceTransformMode == DT_3D) {
    // Calculate 3D Radial transform
	  min = double(EDT_MAX_DISTANCE_SQUARED);
	  for(x=0;x<nx;x++){
		  for(y=0;y<ny;y++){
			  for(z=0;z<nz;z++){
				if(this->_Input->GetAsDouble(x,y,z) < min)
					min = this->_Input->GetAsDouble(x,y,z);
			  }
		  }
	  }
	  for(x=0;x<nx;x++){
		  for(y=0;y<ny;y++){
			  for(z=0;z<nz;z++){
				  this->_Output->PutAsDouble(x,y,z,this->_Input->GetAsDouble(x,y,z)-min);
			  }
		  }
	  }
  } else {
	  // Calculate 2D Radial transform	  
	  for(z=0;z<nz;z++){
		  min = double(EDT_MAX_DISTANCE_SQUARED);
		  for(y=0;y<ny;y++){
			  for(x=0;x<nz;x++){
				if(this->_Input->GetAsDouble(x,y,z) < min)
					min = this->_Input->GetAsDouble(x,y,z);
			  }
		  }
		  for(y=0;y<ny;y++){
			  for(x=0;x<nx;x++){
				  this->_Output->PutAsDouble(x,y,z,this->_Input->GetAsDouble(x,y,z)-min);
			  }
		  }
	  }
  }

  // Do the final cleaning up
  this->Finalize();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void EuclideanDistanceTransform<VoxelType>::TRadial()
{
  int nx, ny, nz, x, y, z;
  double min;

  // Calculate image dimensions
  nx = this->_Input->X();
  ny = this->_Input->Y();
  nz = this->_Input->Z();

  if (this->_distanceTransformMode == DT_3D) {
    // Calculate 3D Radial transform
	  min = double(EDT_MAX_DISTANCE_SQUARED);
	  for(x=0;x<nx;x++){
		  for(y=0;y<ny;y++){
			  for(z=0;z<nz;z++){
				if(this->_Input->GetAsDouble(x,y,z) < min)
					min = this->_Input->GetAsDouble(x,y,z);
			  }
		  }
	  }
	  for(x=0;x<nx;x++){
		  for(y=0;y<ny;y++){
			  for(z=0;z<nz;z++){
				  this->_Output->PutAsDouble(x,y,z,this->_Input->GetAsDouble(x,y,z)-min
					  + abs(this->_Input->GetAsDouble(x,y,z))/2);
			  }
		  }
	  }
  } else {
	  // Calculate 2D Radial transform	  
	  for(z=0;z<nz;z++){
		  min = double(EDT_MAX_DISTANCE_SQUARED);
		  for(y=0;y<ny;y++){
			  for(x=0;x<nz;x++){
				if(this->_Input->GetAsDouble(x,y,z) < min)
					min = this->_Input->GetAsDouble(x,y,z);
			  }
		  }
		  for(y=0;y<ny;y++){
			  for(x=0;x<nx;x++){
				  this->_Output->PutAsDouble(x,y,z,this->_Input->GetAsDouble(x,y,z)-min
					  + abs(this->_Input->GetAsDouble(x,y,z))/2);
			  }
		  }
	  }
  }

  // Do the final cleaning up
  this->Finalize();
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class EuclideanDistanceTransform<float>;
template class EuclideanDistanceTransform<double>;


} // namespace mirtk
