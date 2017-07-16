/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
 * Copyright 2013-2017 Andreas Schuh
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

#ifndef MIRTK_GenericRegistrationFilter_H
#define MIRTK_GenericRegistrationFilter_H

#include "mirtk/RegistrationFilter.h"

#include "mirtk/Math.h"
#include "mirtk/Point.h"
#include "mirtk/Vector3D.h"
#include "mirtk/EventDelegate.h"

#include "mirtk/BaseImage.h"
#include "mirtk/GenericImage.h"
#include "mirtk/RegisteredImage.h"
#include "mirtk/InterpolationMode.h"
#include "mirtk/ExtrapolationMode.h"

#include "mirtk/SimilarityMeasure.h"
#include "mirtk/PointSetDistanceMeasure.h"
#include "mirtk/ConstraintMeasure.h"

#include "mirtk/LocalOptimizer.h"
#include "mirtk/RegistrationEnergy.h"

#include "mirtk/Transformation.h"
#include "mirtk/TransformationModel.h"

#include "mirtk/RegistrationConfig.h"
#if MIRTK_Registration_WITH_PointSet
#  include "vtkSmartPointer.h"
#  include "vtkPointSet.h"
#else
class vtkPointSet;
#endif


namespace mirtk {


class RegisteredPointSet;
class RegistrationEnergyParser;
class GenericRegistrationLogger;
class GenericRegistrationDebugger;
class HomogeneousTransformation;
class FreeFormTransformation;


const int MAX_NO_RESOLUTIONS = 10;


/**
 * Generic registration filter
 */
class GenericRegistrationFilter : public RegistrationFilter
{
  mirtkObjectMacro(GenericRegistrationFilter);

  friend class RegistrationEnergyParser;
  friend class GenericRegistrationLogger;
  friend class GenericRegistrationDebugger;

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Type of resolution pyramid images
  typedef RegisteredImage::InputImageType   ResampledImageType;

  /// List type storing images for one resolution pyramid level
  typedef Array<ResampledImageType>   ResampledImageList;

  /// Scalar type of resolution pyramid images
  typedef ResampledImageType::VoxelType   VoxelType;

  /// Type of cached displacement field
  typedef RegisteredImage::DisplacementImageType   DisplacementImageType;

  /// Structure storing information about transformation instance
  struct TransformationInfo
  {
    double _Exponent;

    TransformationInfo(double e = .0) : _Exponent(e) {}

    bool operator ==(const TransformationInfo &other) const
    {
      return fequal(_Exponent, other._Exponent, 1e-3);
    }

    bool operator !=(const TransformationInfo &other) const
    {
      return !(*this == other);
    }

    static TransformationInfo Identity()
    {
      return TransformationInfo(.0);
    }

    static TransformationInfo Full()
    {
      return TransformationInfo(1.0);
    }

    static TransformationInfo Inverse()
    {
      return TransformationInfo(-1.0);
    }

    bool IsForwardTransformation() const
    {
      return _Exponent > .0;
    }

    bool IsBackwardTransformation() const
    {
      return _Exponent < .0;
    }

    bool IsIdentity() const
    {
      return fequal(_Exponent, .0, 1e-3);
    }

    operator bool() const
    {
      return !IsIdentity();
    }
  };

  /// Structure storing information of a used image similarity term parsed
  /// from the registration energy function string
  struct ImageSimilarityInfo
  {
    SimilarityMeasure  _Measure;              ///< Type of similarity measure
    bool               _DefaultSign;          ///< Whether to use default sign of similarity
    string             _Name;                 ///< Name of similarity term
    double             _Weight;               ///< Weight of similarity term
    int                _TargetIndex;          ///< Index of target image
    TransformationInfo _TargetTransformation; ///< Target transformation identifier
    int                _SourceIndex;          ///< Index of source image
    TransformationInfo _SourceTransformation; ///< Source transformation identifier

    bool IsSymmetric() const
    {
      return _TargetTransformation && _SourceTransformation;
    }
  };

  /// Structure storing information of a used point set distance term parsed
  /// from the registration energy function string
  struct PointSetDistanceInfo
  {
    PointSetDistanceMeasure _Measure;              ///< Measure of polydata distance
    bool                    _DefaultSign;          ///< Whether to use default sign of distance measure
    string                  _Name;                 ///< Name of polydata distance term
    double                  _Weight;               ///< Weight of polydata distance term
    int                     _TargetIndex;          ///< Index of target data set
    TransformationInfo      _TargetTransformation; ///< Target transformation identifier
    int                     _SourceIndex;          ///< Index of source data set
    TransformationInfo      _SourceTransformation; ///< Source transformation identifier

    bool IsSymmetric() const
    {
      return _TargetTransformation && _SourceTransformation;
    }
  };

  /// Structure storing information of a used point set constraint term parsed
  /// from the registration energy function string
  struct PointSetConstraintInfo
  {
    EnergyMeasure      _Measure;          ///< Type of internal forces
    string             _Name;             ///< Name of constraint
    double             _Weight;           ///< Weight of constraint
    int                _PointSetIndex;    ///< Index of input point set object
    int                _RefImageIndex;
    int                _RefPointSetIndex;
    TransformationInfo _Transformation;   ///< Point set transformation identifier
  };

  /// Structure storing information of a used constraint term parsed from the
  /// registration energy function string
  struct ConstraintInfo
  {
    ConstraintMeasure _Measure; ///< Type of constraint
    string            _Name;    ///< Name of constraint
    double            _Weight;  ///< Weight of constraint
  };

  /// Structure storing information about cached displacements
  struct DisplacementInfo
  {
    int                   _DispIndex;      ///< Index of cached displacement field
    double                _InputTime;      ///< Time of untransformed data
    ImageAttributes       _Domain;         ///< Domain on which it is defined
    const Transformation *_Transformation; ///< Corresponding transformation instance
  };

  /// Structure storing information about transformed output point set
  struct PointSetOutputInfo
  {
    int                _InputIndex;
    bool               _InitialUpdate;
    TransformationInfo _Transformation;
  };

  // ---------------------------------------------------------------------------
  // Attributes

  /// Number of resolution levels
  mirtkPublicAttributeMacro(int, NumberOfLevels);

  /// Level at which to stop multi-resolution optimization
  mirtkPublicAttributeMacro(int, FinalLevel);

  /// Multi-level transformation mode
  mirtkPublicAttributeMacro(MFFDMode, MultiLevelMode);

  /// Transformation model
  mirtkPublicAttributeMacro(Array<enum TransformationModel>, TransformationModel);

  /// Default image interpolation mode
  mirtkPublicAttributeMacro(enum InterpolationMode, DefaultInterpolationMode);

  /// Default image extrapolation mode
  mirtkPublicAttributeMacro(enum ExtrapolationMode, DefaultExtrapolationMode);

  /// Rescale input image intensities from [imin, imax] to [0, omax]
  ///
  /// Here [imin, imax] is the intensity range of the respective input image
  /// excluding the background (if specified) and omax is the maximum intensity
  /// of the rescaled intensities. This maximum value influences the relative
  /// weighting of image similarity gradient vs. the gradient of constraint terms
  /// such as the bending energy.
  ///
  /// This normalization ensures that the magnitude of image derivatives used
  /// for energy gradient computation have comparable magnitude in order to
  /// prevent one image to have stronger influence than another. This is more
  /// important for symmetric, inverse consistent, or multi-modal settings.
  ///
  /// By default, this parameter is set to inf and no rescaling is done.
  mirtkPublicAttributeMacro(double, MaxRescaledIntensity);

  /// Whether to precompute image derivatives or compute them on the fly
  mirtkPublicAttributeMacro(bool, PrecomputeDerivatives);

  /// Default similarity measure
  mirtkPublicAttributeMacro(enum SimilarityMeasure, SimilarityMeasure);

  /// Default polydata distance measure
  mirtkPublicAttributeMacro(enum PointSetDistanceMeasure, PointSetDistanceMeasure);

  /// Optimization method
  mirtkPublicAttributeMacro(enum OptimizationMethod, OptimizationMethod);

  /// Normalize weights of energy function terms
  mirtkPublicAttributeMacro(bool, NormalizeWeights);

  /// Initial guess of optimal transformation
  mirtkPublicAggregateMacro(const Transformation, InitialGuess);

  /// Target transformation to be approximated by output transformation
  mirtkPublicAggregateMacro(const Transformation, TargetTransformation);

  /// Name of target transformation MSDE term
  mirtkPublicAttributeMacro(string, TargetTransformationErrorName);

  /// Weight of target transformation MSDE term
  mirtkPublicAttributeMacro(double, TargetTransformationErrorWeight);

  /// Whether to merge initial global transformation into (local) transformation
  mirtkPublicAttributeMacro(bool, MergeGlobalAndLocalTransformation);

  /// Whether to allow x coordinate transformation
  mirtkPublicAttributeMacro(bool, RegisterX);

  /// Whether to allow y coordinate transformation
  mirtkPublicAttributeMacro(bool, RegisterY);

  /// Whether to allow z coordinate transformation
  mirtkPublicAttributeMacro(bool, RegisterZ);

  /// Enforce Dirichlet boundary condition on FFD transformations
  ///
  /// When this option is enabled, the status of control points at the
  /// boundary of the finite FFD lattice is set to Passive and the parameters
  /// of these control points set to zero. This is always the case for FFDs
  /// whose parameters are control point displacements. In case of FFDs
  /// parameterized by (stationary) velocities, the default extrapolation
  /// mode is nearest neighbor, however, and a layer of passive CPs with
  /// constant value is needed if the velocity should be forced to zero
  /// outside the finite domain on which the velocity field is defined.
  /// Alternatively, FFD extrapolation mode "Const" can be used.
  mirtkPublicAttributeMacro(bool, DirichletBoundaryCondition);

  /// Mask which defines where to evaluate the energy function
  mirtkPublicAggregateMacro(BinaryImage, Domain);

  /// Whether to adaptively remesh surfaces before each gradient step
  mirtkPublicAttributeMacro(bool, AdaptiveRemeshing);

protected:

  /// Common attributes of (untransformed) input target data sets
  ///
  /// These attributes are in particular used to initialize the control point
  /// grid of a free-form deformation such that the grid is large enough to
  /// be valid for every target point for which the transformation will be
  /// evaluated. It is therefore computed from the attributes of the input
  /// data sets rather than the downsampled data. If a _Domain mask is set,
  /// the attributes of this mask are copied instead.
  struct ImageAttributes _RegistrationDomain;

  Array<Point>               _Centroid;                       ///< Centroids of images
  Point                      _TargetOffset;                   ///< Target origin offset
  Point                      _SourceOffset;                   ///< Source origin offset
  Array<const BaseImage *>   _Input;                          ///< Input images
  Array<ResampledImageList>  _Image;                          ///< Resolution pyramid
  Array<BinaryImage *>       _Mask;                           ///< Domain on which to evaluate similarity
  Transformation            *_Transformation;                 ///< Current estimate
  Array<TransformationInfo>  _TransformationInfo;             ///< Meta-data of partial transformation
  Array<Transformation *>    _TransformationInstance;         ///< Partial transformations
  Array<enum InterpolationMode> _InterpolationMode;           ///< Interpolation mode of each input image
  Array<enum ExtrapolationMode> _ExtrapolationMode;           ///< Extrapolation mode of each input image
  RegistrationEnergy         _Energy;                         ///< Registration energy
  LocalOptimizer            *_Optimizer;                      ///< Used optimizer
  enum TransformationModel   _CurrentModel;                   ///< Current transformation model
  int                        _CurrentLevel;                   ///< Current resolution level
  EventDelegate              _EventDelegate;                  ///< Forwards optimization events to observers
  string                     _EnergyFormula;                  ///< Registration energy formula as string
  Array<ImageSimilarityInfo> _ImageSimilarityInfo;            ///< Parsed similarity measure(s)
  Array<ConstraintInfo>      _ConstraintInfo;                 ///< Parsed constraint(s)
  Array<Vector3D<double> >   _Resolution[MAX_NO_RESOLUTIONS]; ///< Image resolution in mm
  Array<double>              _MinEdgeLength[MAX_NO_RESOLUTIONS]; ///< Minimum edge length in mm
  Array<double>              _MaxEdgeLength[MAX_NO_RESOLUTIONS]; ///< Maximum edge length in mm
  int                        _UseGaussianResolutionPyramid;   ///< Whether resolution levels correspond to a Gaussian pyramid
  Array<double>              _Blurring[MAX_NO_RESOLUTIONS];   ///< Image blurring value
  double                     _DefaultBackground;              ///< Default background value
  Array<double>              _Background;                     ///< Image background value
  bool                       _DownsampleWithPadding;          ///< Whether to take background into account
                                                              ///< during initialization of the image pyramid
  bool                       _CropPadImages;                  ///< Whether to crop/pad input images
  int                        _CropPadFFD;                     ///< Whether to crop/pad FFD lattice

#if MIRTK_Registration_WITH_PointSet
  Array<double>                   _PointSetTime;              ///< Time point of input points, curves, and/or surfaces
  Array<vtkSmartPointer<vtkPointSet> > _PointSetInput;        ///< Input points, curves, and/or surfaces
  Array<Array<vtkSmartPointer<vtkPointSet> > > _PointSet;     ///< Remeshed/-sampled points, curves, and/or surfaces
  Array<RegisteredPointSet *>     _PointSetOutput;            ///< Output points, curves, and/or surfaces
  Array<PointSetOutputInfo>       _PointSetOutputInfo;        ///< Meta-data of output points, curves, and/or surfaces
  Array<PointSetDistanceInfo>     _PointSetDistanceInfo;      ///< Parsed point set distance measure(s)
  Array<PointSetConstraintInfo>   _PointSetConstraintInfo;    ///< Parsed point set constraint measure(s)
#else // MIRTK_Registration_WITH_PointSet
  Array<double>                   _PointSetTime;              ///< Unused
  Array<void *>                   _PointSetInput;             ///< Unused
  Array<Array<void *> >           _PointSet;                  ///< Unused
  Array<void *>                   _PointSetOutput;            ///< Unused
  Array<PointSetOutputInfo>       _PointSetOutputInfo;        ///< Unused
  Array<PointSetDistanceInfo>     _PointSetDistanceInfo;      ///< Unused
  Array<PointSetConstraintInfo>   _PointSetConstraintInfo;    ///< Unused
#endif // MIRTK_Registration_WITH_PointSet

  int    _Centering             [MAX_NO_RESOLUTIONS];      ///< Whether to center foreground (if applicable)
  double _MinControlPointSpacing[MAX_NO_RESOLUTIONS][4];   ///< Control point spacing for FFDs
  double _MaxControlPointSpacing[MAX_NO_RESOLUTIONS][4];   ///< Control point spacing for FFDs
  bool   _Subdivide             [MAX_NO_RESOLUTIONS][4];   ///< Whether to subdivide FFD

  /// Parameters not considered by the registration filter itself
  ///
  /// These parameters are passed on to the sub-modules of the registration
  /// filter such as the image similarity measure(s), the registration
  /// constraint(s), and the optimizer. This way, the registration filter does
  /// not need to know which parameters its sub-modules accept and is better
  /// decoupled from the implementation of the respective image similarities,
  /// constraints, and optimizers.
  ParameterList _Parameter[MAX_NO_RESOLUTIONS];

  /// Delegate of pre-update function
  RegistrationEnergy::PreUpdateFunctionType _PreUpdateDelegate;

  /// Info of cached displacement fields
  Array<DisplacementInfo>        _DisplacementInfo;
  Array<DisplacementImageType *> _DisplacementField;

  // ---------------------------------------------------------------------------
  // Access to managed data objects and their attributes

protected:

  /// Determine number of temporal frames and temporal sampling attributes
  virtual int NumberOfFrames(double * = NULL, double * = NULL, double * = NULL) const;

  /// Get attributes of n-th input image at specified resolution level
  virtual struct ImageAttributes ImageAttributes(int, int = -1) const;

  /// Common attributes of (untransformed) input target data sets at specified resolution level
  virtual struct ImageAttributes RegistrationDomain(int = -1) const;

  /// Average resolution of (untransformed) target data sets at given level
  virtual Vector3D<double> AverageOutputResolution(int = -1) const;

  /// Get (partial/inverse) output transformation
  virtual Transformation *OutputTransformation(TransformationInfo);

  /// Initialize output image corresponding to registered (transformed) image
  virtual void SetInputOf(RegisteredImage *, const struct ImageAttributes &, int, TransformationInfo);

  /// Get (new) registered output point set
  ///
  /// \note Use only when MIRTK_Registration_WITH_PointSet is 1.
  RegisteredPointSet *OutputPointSet(int, double, TransformationInfo);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
private:

  /// Copy constructor
  /// \note Intentionally not implemented
  GenericRegistrationFilter(const GenericRegistrationFilter &);

  /// Assignment operator
  /// \note Intentionally not implemented
  void operator =(const GenericRegistrationFilter &);

public:

  /// Constructor
  GenericRegistrationFilter();

  /// Reset filter settings, but keep input
  virtual void Reset();

  /// Reset filter settings and input
  virtual void Clear();

  /// Destructor
  virtual ~GenericRegistrationFilter();

  // ---------------------------------------------------------------------------
  // Input images

  /// Set input images of the registration filter
  void Input(const BaseImage *, const BaseImage *);

  /// Set input images of the registration filter
  void Input(int, const BaseImage **);

  /// Set input images of the registration filter
  template <class TVoxel> void Input(int, const GenericImage<TVoxel> **);

  /// Add filter input image
  ///
  /// For registration filters which support multiple target/source images
  /// such as multiple channels and/or frames of a temporal image sequence.
  /// How the transformation is applied to which inputs depends on the
  /// respective transformation model and image registration method.
  void AddInput(const BaseImage *);

  /// Number of input images
  int NumberOfImages() const;

  /// Number of required input images
  /// \note Only valid after ParseEnergyFormula has been called!
  int NumberOfRequiredImages() const;

  /// Determine whether the specified input image either remains untransformed,
  /// or is transformed by the inverse transformation or a part of it such
  /// as in case of an inverse consistent and/or symmetric energy function.
  bool IsTargetImage(int) const;

  /// Determine whether the specified input image will be transformed by the
  /// forward transformation or part of it such as in case of a symmetric energy.
  bool IsSourceImage(int) const;

  /// Determine whether the specified input image remains untransformed.
  bool IsFixedImage(int) const;

  /// Determine whether the specified input image will be transformed.
  bool IsMovingImage(int) const;

  /// Set common/default interpolation mode for all input images
  void InterpolationMode(enum InterpolationMode);

  /// Set interpolation mode of n-th input image
  void InterpolationMode(int, enum InterpolationMode);

  /// Get interpolation mode of n-th input image or default mode, respectively
  enum InterpolationMode InterpolationMode(int = -1) const;

  /// Set common/default extrapolation mode for all input images
  void ExtrapolationMode(enum ExtrapolationMode);

  /// Set extrapolation mode of n-th input image
  void ExtrapolationMode(int, enum ExtrapolationMode);

  /// Get extrapolation mode of n-th input image or default mode, respectively
  enum ExtrapolationMode ExtrapolationMode(int = -1) const;

  /// Get background value of n-th input image (after registration done)
  double BackgroundValue(int) const;

  // ---------------------------------------------------------------------------
  // Input simplicial complexes (points, curves, surfaces, tetrahedral meshes)

  /// Set input point sets of the registration filter
  ///
  /// \note Use only when MIRTK_Registration_WITH_PointSet is 1.
  void Input(vtkPointSet *, vtkPointSet *, double = .0, double = 1.0);

  /// Set input point set of the registration filter
  ///
  /// \note Use only when MIRTK_Registration_WITH_PointSet is 1.
  void Input(int, vtkPointSet **, double * = NULL);

  /// Add filter input point set
  ///
  /// \note Use only when MIRTK_Registration_WITH_PointSet is 1.
  void AddInput(vtkPointSet *, double = .0);

  /// Number of input point sets
  int NumberOfPointSets() const;

  /// Number of required input point sets
  /// \note Only valid after ParseEnergyFormula has been called!
  int NumberOfRequiredPointSets() const;

  /// Determine whether the specified input point set will be transformed by the
  /// forward transformation or part of it such as in case of a symmetric energy.
  bool IsTargetPointSet(int) const;

  /// Determine whether the specified input point set either remains untransformed,
  /// or is transformed by the inverse transformation or a part of it such
  /// as in case of an inverse consistent and/or symmetric energy function.
  bool IsSourcePointSet(int) const;

  /// Determine whether the specified input point set remains untransformed.
  bool IsFixedPointSet(int) const;

  /// Determine whether the specified input point set will be transformed.
  bool IsMovingPointSet(int) const;

  // ---------------------------------------------------------------------------
  // Parameter

  using RegistrationFilter::Read;
  using RegistrationFilter::Parameter;

  /// Set (single) transformation model
  virtual void TransformationModel(enum TransformationModel);

  /// Parse registration energy function
  /// \note Optional to call explicitly as this method is called by GuessParameter.
  virtual void ParseEnergyFormula(int = -1, int = -1, int = -1);

  /// Guess proper setting for any yet unset parameter
  /// \note Optional to call explicitly as this method is called by Run.
  virtual void GuessParameter();

  /// Read registration parameters from input stream
  virtual bool Read(istream &, bool = false);

  /// Set named parameter from value as string
  virtual bool Set(const char *, const char *);

  /// Set named parameter from value as string
  virtual bool Set(const char *, const char *, int);

  /// Get parameters as key/value as string map
  virtual ParameterList Parameter() const;

  /// Get parameters as key/value as string map
  virtual ParameterList Parameter(int) const;

  /// Write registration parameters to file
  virtual void Write(const char *) const;

  // ---------------------------------------------------------------------------
  // Execution

  /// Make an initial guess of the (global) output transformation
  virtual Transformation *MakeInitialGuess();

  /// Run the multi-level registration
  virtual void Run();

  // ---------------------------------------------------------------------------
  // Implementation - the following functions can be overridden in subclasses

protected:

  /// Whether current level is the initial resolution level
  bool AtInitialLevel() const;

  /// Whether current level is the final resolution level
  bool AtFinalLevel() const;

  /// Run multi-resolution registration
  virtual void MultiResolutionOptimization();

  /// Initialize registration at current resolution
  virtual void Initialize();

  /// Initialize image resolution pyramid
  virtual void InitializePyramid();

  /// Remesh/-sample input point sets
  virtual void InitializePointSets();

  /// Type of output transformation of sub-registration at current resolution
  virtual enum TransformationType TransformationType();

  /// Initialize new transformation instance
  virtual void InitializeTransformation();

  /// Initialize transformation parameters using provided initial guess
  virtual void ApplyInitialGuess();

  /// Initialize status of linear parameters
  virtual void InitializeStatus(HomogeneousTransformation *);

  /// Initialize status of FFD parameters
  virtual void InitializeStatus(FreeFormTransformation *);

  /// Initialize status of transformation parameters
  virtual void InitializeStatus();

  /// Initialize transformation for sub-registration at current resolution
  virtual void InitializeOutput();

  /// Instantiate new image similarity term for energy function
  /// \note The individual energy terms are destroyed by the energy function!
  virtual void AddImageSimilarityTerm();

  /// Instantiate new point set distance term for energy function
  /// \note The individual energy terms are destroyed by the energy function!
  virtual void AddPointSetDistanceTerm();

  /// Instantiate new point set constraint term for energy function
  /// \note The individual energy terms are destroyed by the energy function!
  virtual void AddPointSetConstraintTerm();

  /// Instantiate new regularization term for energy function
  /// \note The individual energy terms are destroyed by the energy function!
  virtual void AddPenaltyTerm();

  /// Initialize registration energy of registration at current resolution
  virtual void InitializeEnergy();

  /// Initialize optimizer used to solve registration problem
  virtual void InitializeOptimizer();

  /// Finalize registration at current resolution
  virtual void Finalize();

  /// Callback function called by _Energy->Update(bool)
  void PreUpdateCallback(bool);

};

////////////////////////////////////////////////////////////////////////////////
// Inline/template definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <class TVoxel>
void GenericRegistrationFilter::Input(int num, const GenericImage<TVoxel> **image)
{
  _Input.clear();
  for (int n = 0; n < num; ++n) AddInput(image[n]);
}


} // namespace mirtk

#endif // MIRTK_GenericRegistrationFilter_H
