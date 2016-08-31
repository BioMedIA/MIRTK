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

#ifndef MIRTK_MeshSmoothing_H
#define MIRTK_MeshSmoothing_H

#include "mirtk/MeshFilter.h"

#include "mirtk/Array.h"

class vtkDataArray;


namespace mirtk {


/**
 * Smooth scalars and/or points of triangulated surface mesh
 */
class MeshSmoothing : public MeshFilter
{
  mirtkObjectMacro(MeshSmoothing);

  // ---------------------------------------------------------------------------
  // Types

public:

  /// Enumeration of smoothing kernel functions
  enum WeightFunction
  {
    Default,             ///< Undefined weighting function, use default
    Combinatorial,       ///< Uniform node weights / "umbrella operator"
    InverseDistance,     ///< Inverse node distance
    Gaussian,            ///< Gaussian node weights
    AnisotropicGaussian, ///< Anisotropic Gaussian node weights
    NormalDeviation      ///< Weight by cosine of angle made up by normals
  };

  /// List of point data array names
  typedef Array<string> ArrayNames;

  /// Vector of point data arrays to be smoothed
  typedef Array<vtkSmartPointer<vtkDataArray> > DataArrays;

  // ---------------------------------------------------------------------------
  // Attributes

private:

  /// Input point mask, only points with mask value != 0 are modified
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataArray>, Mask);

  /// Number of smoothing iterations
  mirtkPublicAttributeMacro(int, NumberOfIterations);

  /// Relaxation / scale factor for odd iterations
  mirtkPublicAttributeMacro(double, Lambda);

  /// Relaxation / scale factor for even iterations
  ///
  /// By default, when _Mu is set to NaN, a _Mu value equal to _Lambda is used.
  /// A negative _Mu value, greater in magnitude than _Lambda, can be used
  /// to perform a low-pass filtering of the surface mesh as described in:
  /// Taubin, G., Curve and Surface Smoothing without Shrinkage, ICCV 1995.
  mirtkPublicAttributeMacro(double, Mu);

  /// Smoothing kernel parameter
  ///
  /// In case of a Gaussian smoothing kernel, if the sigma value is negative,
  /// the standard deviation of the Gaussian kernel is set to the average
  /// edge length times the absolute value of \c _Sigma. If _Sigma is zero,
  /// the average edge length is used as standard deviation. The _Sigma attribute
  /// is set to the actual used standard deviation after the filter execution.
  ///
  /// In case of the inverse distance weighting, the _Sigma value is added to
  /// the edge length before computing the inverse value.
  mirtkPublicAttributeMacro(double, Sigma);

  /// Smoothing kernel parameter in direction of maximum curvature
  ///
  /// \note The direction of maximum principle curvature is orthogonal to the
  ///       direction in which the surface is most bended! It is the direction
  ///       with the most variance, i.e., along ridges, not orthogonal to these.
  ///
  /// This parameter is only used by the anisotropic Gaussian kernel.
  /// By default, i.e., when _Sigma2 = 0, the standard deviation along the
  /// direction of maximum change is half the standard deviation along the
  /// direction of minimum change. Hence, the surface points or data values are
  /// smoothed less in the direction of maximum change (i.e., maximum curvature).
  /// If the sigma value is negative, the standard deviation is set to the
  /// average edge length times the absolute value of \c _Sigma2.
  ///
  /// When an array of local geometry tensors is used instead of the direction
  /// of minimum and/or maximum change, the default is to use an isotropic
  /// Gaussian kernel in the local coordinate system defined by the tensor.
  /// In this case the axes of the local coordinate system are expected to be
  /// scaled anisotropically as in case of the curvature tensor, for example.
  mirtkPublicAttributeMacro(double, MaximumDirectionSigma);

  /// Smoothing kernel function
  mirtkPublicAttributeMacro(WeightFunction, Weighting);

  /// Name of input point data array with local geometry tensor used for anisotropic smoothing
  ///
  /// For example, the local curvature tensors can be computed using
  /// PolyDataCurvature and used for anisotropic Gaussian smoothing.
  /// The input point data array is then named PolyDataCurvature::TENSOR.
  mirtkPublicAttributeMacro(string, GeometryTensorName);

  /// Name of input point data array with direction along which to smooth less
  /// \note This array is only used if no _GeometryTensorName is specified.
  mirtkPublicAttributeMacro(string, MinimumDirectionName);

  /// Name of input point data array with direction along which to smooth more
  /// \note This array is only used if no _GeometryTensorName is specified.
  mirtkPublicAttributeMacro(string, MaximumDirectionName);

  /// Whether to average values of adjacent nodes only or to
  /// also include the node's values themselves in the average
  ///
  /// \note In case of an InverseDistance node weighting, the values of the
  ///       node itself are only included in the average if _Sigma > .0.
  mirtkPublicAttributeMacro(bool, AdjacentValuesOnly);

  /// Whether to smooth the node positions, i.e., input geometry
  mirtkPublicAttributeMacro(bool, SmoothPoints);

  /// Smooth magnitude of (vector-valued) data arrays
  mirtkPublicAttributeMacro(bool, SmoothMagnitude);

  /// Smooth (vector-valued) data arrays by only considering those adjacent
  /// vectors which have a positive dot product with the data vector at the
  /// current node position
  mirtkPublicAttributeMacro(bool, SignedSmoothing);

  /// Names of input point data arrays to be smoothed
  mirtkPublicAttributeMacro(ArrayNames, SmoothArrays);

  /// Array vtkDataSetAttributes::AttributeTypes
  mirtkAttributeMacro(Array<int>, AttributeTypes);

  /// Input point data arrays to be smoothed
  mirtkAttributeMacro(DataArrays, InputArrays);

  /// Output point data arrays
  mirtkAttributeMacro(DataArrays, OutputArrays);

  /// Verbosity of output messages
  mirtkPublicAttributeMacro(int, Verbose);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const MeshSmoothing &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  MeshSmoothing();

  /// Copy constructor
  MeshSmoothing(const MeshSmoothing &);

  /// Assignment operator
  MeshSmoothing &operator =(const MeshSmoothing &);

  /// Destructor
  virtual ~MeshSmoothing();

  /// Add named point data array to list of arrays to be smoothed
  ///
  /// \param[in] name Name of data array to smooth.
  /// \param[in] attr Array attribute type, see vtkDataSetAttributes::AttributeTypes.
  ///                 When vtkDataSetAttributes::VECTORS, an array with 3 components
  ///                 is treated as 3D direction vectors and the smoothing is
  ///                 performed independent of the sign of the direction.
  void SmoothArray(const char *name, int attr = -1);

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Execute filter
  virtual void Execute();

  /// Finalize filter execution
  virtual void Finalize();

  // ---------------------------------------------------------------------------
  // Alternative VTK-like interface

public:

  /// Set number of smoothing iterations
  mirtkSetMacro(NumberOfIterations, int);

  /// Get number of smoothing iterations
  mirtkGetMacro(NumberOfIterations, int);

  /// Set relaxation factor
  mirtkSetMacro(Lambda, double);

  /// Get relaxation factor
  mirtkGetMacro(Lambda, double);

  /// Set smoothing kernel standard deviation
  mirtkSetMacro(Sigma, double);

  /// Get smoothing kernel standard deviation
  mirtkGetMacro(Sigma, double);

  /// Set smoothing kernel standard deviation in direction of maximum curvature
  mirtkSetMacro(MaximumDirectionSigma, double);

  /// Get smoothing kernel standard deviation in direction of maximum curvature
  mirtkGetMacro(MaximumDirectionSigma, double);

  /// Enable/disable averaging of adjacent node values only
  mirtkOnOffMacro(AdjacentValuesOnly);

  /// Enable/disable smoothing of node positions
  mirtkOnOffMacro(SmoothPoints);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <>
inline bool FromString(const char *str, MeshSmoothing::WeightFunction &value)
{
  const string lstr = ToLower(Trim(str));
  if      (lstr == "combinatorial")       value = MeshSmoothing::Combinatorial;
  else if (lstr == "inversedistance")     value = MeshSmoothing::InverseDistance;
  else if (lstr == "gaussian")            value = MeshSmoothing::Gaussian;
  else if (lstr == "anisotropicgaussian") value = MeshSmoothing::AnisotropicGaussian;
  else {
    value = MeshSmoothing::Default;
    return (lstr == "default");
  }
  return true;
}

// -----------------------------------------------------------------------------
inline void MeshSmoothing::SmoothArray(const char *name, int attr)
{
  _SmoothArrays.push_back(name);
  _AttributeTypes.push_back(attr);
}


} // namespace mirtk

#endif // MIRTK_MeshSmoothing_H