/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkDummyImageSimilarity.h>


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace irtkDummyImageSimilarityUtils {


// -----------------------------------------------------------------------------
// Types
typedef irtkDummyImageSimilarity::VoxelType    VoxelType;
typedef irtkDummyImageSimilarity::GradientType GradientType;

// -----------------------------------------------------------------------------
/// Parallel evaluation of image dissimilarity
struct Evaluate : public irtkVoxelReduction
{
  irtkDummyImageDissimilarity *_Sim;
  double                       _Sum;
  int                          _Cnt;

  Evaluate(irtkDummyImageSimilarity *sim)
  :
    _Sim(sim), _Sum(.0), _Cnt(0)
  {}

  Evaluate(const Evaluate &rhs)
  :
    _Sim(rhs._Sim), _Sum(rhs._Sum), _Cnt(rhs._Cnt)
  {}

  void split(Evaluate &lhs)
  {
    _Sum = .0;
    _Cnt =  0;
  }

  void join(Evaluate &rhs)
  {
    _Sum += rhs._Sum;
    _Cnt += rhs._Cnt;
  }

  void operator ()(int i, int j, int k, int, VoxelType *t, VoxelType *s)
  {
    if (_Sim->IsForeground(i, j, k)) {
      // Add dissimilarity of intensity values *t and *s
      _Sum += static_cast<double>(pow(*t - *s, 2));
      // Increase counter of voxels within foreground/overlap region
      ++_Cnt;
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate voxel-based gradient of image dissimilarity measure
struct EvaluateGradient : public irtkVoxelFunction
{
  irtkDummyImageSimilarity *_Sim;

  EvaluateGradient(irtkDummyImageSimilarity *sim)
  :
    _Sim(sim)
  {}

  void operator ()(int i, int j, int k, int, const VoxelType *t, const VoxelType *s, GradientType *g)
  {
    if (_Sim->IsForeground(i, j, k)) {
      // Evaluate gradient of dissimilarity measure for given pair of image
      // intensity values at voxel with indices (i, j, k)
      *g = -2.0 * static_cast<double>(*t - *s);
    } else {
      // Set gradient to zero outside of foreground/overlap region
      *g = .0;
    }
  }
};


} // namespace irtkDummyImageSimilarityUtils
using namespace irtkDummyImageSimilarityUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkDummyImageSimilarity
::irtkDummyImageSimilarity(const char *name)
:
  irtkImageSimilarity(name),
  _Option(false), _Parameter(.0)
{
}

// -----------------------------------------------------------------------------
irtkDummyImageSimilarity
::irtkDummyImageSimilarity(const irtkDummyImageSimilarity &other)
:
  irtkImageSimilarity(other),
  _Option(other._Option), _Parameter(other._Parameter)
{
}

// -----------------------------------------------------------------------------
irtkDummyImageSimilarity::~irtkDummyImageSimilarity()
{
}

// =============================================================================
// Parameter
// =============================================================================

// -----------------------------------------------------------------------------
irtkDummyImageSimilarity::Set(const char *param, const char *value)
{
  // The irtkEnergyTerm::ParameterNameWithoutPrefix function strips of any
  // common parameter prefix. For example, if param is "Image dissimilarity weight",
  // the returned string is "Weight". The prefix "Image dissimilarity " is added
  // by the superclass irtkImageSimilarity to the list of recognized prefixes.
  // The name of the image similarity term is always considered as parameter
  // name prefix such that the line "<dissimilarity> <parameter> = <value>"
  // in the configuration file identifies this particular instance named
  // "<dissimilarity>". The parameter name is commonly capital with otherwise
  // only lowercase letters, numbers, spaces, and special characters.
  string name = ParameterNameWithoutPrefix(param);
 
  // Set the double parameter
  if (name == "Parameter") {
    // Return true if conversion of string to double was successful.
    // Further, require in this example that the parameter value be positive.
    return FromString(value, _Parameter) && _Parameter > .0;
  }
  if (name == "Option") {
    // Return true if conversion of string to bool was successful.
    // Valid value strings are:
    // - true:  "1" "true"  "True"  "yes" "Yes" "on"  "On"
    // - false: "0" "false" "False" "no"  "No"  "off" "Off"
    return FromString(value, _Option);
  }

  // Not a parameter of this specialized class, but maybe one of the base class
  return irtkImageSimilarity::Set(param, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkDummyImageSimilarity::Parameter() const
{
  // Get list of parameter/value pairs from base class
  irtkParameterList params = irtkImageSimilarity::Parameter();
  // Insert own additional parameter/value pairs if instance has a name
  if (!_Name.empty()) {
    Insert(params, _Name + " parameter", ToString(_Parameter));
    Insert(params, _Name + " option",    ToString(_Option));
  }
  return params;
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void irtkDummyImageSimilarity::Initialize()
{
  // Initialize base class
  irtkImageSimilarity::Initialize();

  // Check that parameters are set properly, possibly choose some default
  // parameter values otherwise or print an error and exit with code 1.

  // Do further necessary initialization after input was set. Note that the
  // input images are *not* yet updated at this point. Only _Target->InputImage()
  // and _Source->InputImage() are valid.
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkDummyImageSimilarity::Update(bool gradient)
{
  // Upate base class and moving image(s)
  irtkImageSimilarity::Update(gradient);

  // Update internal state if necessary because either one of the input images
  // may have changed. This function is called before the gradient evalution
  // of the registration energy function with gradient == true. During the
  // line search with fixed energy function gradient, the parameter
  // gradient == false as only the new image dissimilarity of fixed and
  // moving images has to be re-evaluated to check for improvement.
}

// -----------------------------------------------------------------------------
void irtkDummyImageSimilarity::Exclude(const blocked_range3d<int> &region)
{
  // Not implemented in this example, see irtkSumOfSquaredIntensityDifferences.
  // Instead, the irtkDummyImageSimilarity::Evaluate function always evalutes
  // the dissimilarity of all pairwise voxel intensities.
}

// -----------------------------------------------------------------------------
void irtkDummyImageSimilarity::Include(const blocked_range3d<int> &region)
{
  // Not implemented in this example, see irtkSumOfSquaredIntensityDifferences.
  // Instead, the irtkDummyImageSimilarity::Evaluate function always evalutes
  // the dissimilarity of all pairwise voxel intensities.
}

// -----------------------------------------------------------------------------
double irtkDummyImageSimilarity::Evaluate()
{
  // Parallel evaluation of image dissimilarity (see irtkVoxelFunction.h)
  // Alternatively, simply implement the serial evaluation of the dissimilarity
  // measure between image _Target and _Source. Both are resampled on the image
  // grid defined by the irtkImageAttriutes _Domain, i.e.,
  // - _Target->Attributes() == _Domain
  // - _Source->Attributes() == _Domain
  // As the registration is by default formulated as minimization problem,
  // this function should return a high value for dissimilar images
  // and a low value for similar images, where the lowest value (e.g. zero)
  // is returned for identical images.
  irtkDummyImageSimilarityUtils::Evaluate eval(this);
  ParallelForEachVoxel(_Domain, _Target, _Source, eval);
  return eval._Cnt ? eval._Sum / eval._Cnt : .0;
}

// -----------------------------------------------------------------------------
bool irtkDummyImageSimilarity
::NonParametricGradient(const irtkRegisteredImage *, GradientImageType *)
{
  // No analytic gradient computation implemented by this image dissimilarity.
  // Therefore, false is returned to inform the base class function
  // irtkImageSimilarity::EvaluateGradient that it should approximate the
  // gradient using a finite difference method which only requires the
  // function irtkDummyImageSimilarity::Evaluate to be implemented.
  // For a more efficient approximate gradient computation, implement also
  // irtkDummyImageSimilarity::Exclude and irtkDummyImageSimilarity::Include.
  // This requires to save some internal state of the image dissimilarity to
  // enable a partial evaluation only. See irtkSumOfSquaredIntensityDifferences
  // for a simple example.
  return false;
}
