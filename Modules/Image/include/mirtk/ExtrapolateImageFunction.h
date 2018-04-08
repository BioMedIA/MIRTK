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

#ifndef MIRTK_ExtrapolateImageFunction_H
#define MIRTK_ExtrapolateImageFunction_H

#include "mirtk/ExtrapolationMode.h"
#include "mirtk/Vector.h"
#include "mirtk/VoxelCast.h"
#include "mirtk/BaseImage.h"
#include "mirtk/ImageFunction.h"


namespace mirtk {


/**
 * Abstract base class for discrete image extrapolation
 *
 * This abstract base class defines a common interface for image functions
 * which provide a value at each discrete image location, especially also
 * those which are outside the finite domain on which the image is defined.
 * At these outside voxels, the image function extrapolates the image values,
 * thus extending an image function to an infinite discrete lattice. An image
 * extrapolate image function is mainly used indirectly through an instance of
 * InterpolateImageFunction. Image values at continuous voxel positions
 * are interpolated by the interpolate image function using the extrapolated
 * image values at discrete lattice points.
 */
class ExtrapolateImageFunction : public ImageFunction
{
  mirtkAbstractMacro(ExtrapolateImageFunction);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  ExtrapolateImageFunction();

public:

  /// Destructor
  virtual ~ExtrapolateImageFunction();

  /// Construct extrapolator or return NULL if \p mode is \c Extrapolation_None
  static ExtrapolateImageFunction *New(ExtrapolationMode mode,
                                       const BaseImage * = NULL);

  // ---------------------------------------------------------------------------
  // Initialization

  /// Get extrapolation mode corresponding to this extrapolator
  virtual enum ExtrapolationMode ExtrapolationMode() const = 0;

  // ---------------------------------------------------------------------------
  // Emulation of selected BaseImage functions

  // The following wrapper functions enable the subclasses of
  // InterpolateImageFunction to use a template function for the
  // implementation of EvaluateInside and EvaluateOutside where the only
  // difference is the image type used to access the image values.

  const ImageAttributes &Attributes() const; ///< Attributes of input image

  int X() const; ///< Size of input image in x dimension
  int Y() const; ///< Size of input image in y dimension
  int Z() const; ///< Size of input image in z dimension
  int T() const; ///< Size of input image in t dimension
  int N() const; ///< Number of vector components per image voxel

  double XSize() const; ///< Voxel size of input image in x dimension
  double YSize() const; ///< Voxel size of input image in y dimension
  double ZSize() const; ///< Voxel size of input image in z dimension
  double TSize() const; ///< Voxel size of input image in t dimension

  /// Whether the extrapolated value is identical to the background
  /// value of the input image. If the input uses a background mask
  /// instead, only subclasses of IndexExtrapolateImageFunction
  /// can determine whether the voxel at the transformed location
  /// is in the foreground region.
  virtual bool IsForeground(int, int, int = 0, int = 0) const;

  // ---------------------------------------------------------------------------
  // Image function interface

  /// Get scalar image value at the nearest discrete image location
  double Evaluate(double, double, double = 0, double = 0) const;

  /// Get scalar image value at the nearest discrete image location
  virtual double Evaluate(double, double, double = 0, double = 0);

  // ---------------------------------------------------------------------------
  // Extrapolation of discrete image function

  /// Get scalar image value at an arbitrary discrete image location
  double GetAsDouble(int) const;

  /// Get scalar image value at an arbitrary discrete image location
  virtual double GetAsDouble(int, int, int = 0, int = 0) const = 0;

  /// Get vector image value at an arbitrary discrete image location
  Vector GetAsVector(int) const;

  /// Get vector image value at an arbitrary discrete image location
  Vector GetAsVector(int, int, int = 0, int = 0) const;

  /// Get vector image value at an arbitrary discrete image location
  void GetAsVector(Vector &, int) const;

  /// Get vector image value at an arbitrary discrete image location
  virtual void GetAsVector(Vector &, int, int, int = 0, int = 0) const = 0;

};

////////////////////////////////////////////////////////////////////////////////
// Generic extrapolation interface
////////////////////////////////////////////////////////////////////////////////

/**
 * Generic base class for image extrapolation functions
 *
 * This base class is templated over the type of the input image to be extrapolated
 * and thus subclasses can make use of image voxel type specific (non-virtual)
 * getters to access the image data. No conversion of the voxel type to a uniform
 * voxel type such as double scalar or vector type is required when the type of
 * the image to be extrapolated is known. Otherwise, if the image type is unknown,
 * use the abstract ExtrapolateImageFunction interface instead.
 *
 * \sa ExtrapolateImageFunction
 */
template <class TImage>
class GenericExtrapolateImageFunction : public ExtrapolateImageFunction
{
  mirtkAbstractMacro(GenericExtrapolateImageFunction);

  // ---------------------------------------------------------------------------
  // Types

public:

  typedef TImage                        ImageType; ///< Input image type
  typedef typename ImageType::VoxelType VoxelType; ///< Input voxel type
  typedef typename ImageType::RealType  RealType;  ///< Compatible floating-point type

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  GenericExtrapolateImageFunction();

public:

  /// Destructor
  virtual ~GenericExtrapolateImageFunction();

  /// Construct extrapolator or return NULL if \p mode is \c Extrapolation_None
  static GenericExtrapolateImageFunction *New(enum ExtrapolationMode,
                                              const TImage * = NULL);

  // ---------------------------------------------------------------------------
  // Initialization

  /// Set input image
  virtual void Input(const BaseImage *);

  /// Get input image
  const ImageType *Input() const;

  /// Initialize image function
  virtual void Initialize();

  // ---------------------------------------------------------------------------
  // Extrapolation of discrete image function

  /// Get image value at an arbitrary discrete image location
  VoxelType Get(int) const;

  /// Get image value at an arbitrary discrete image location
  virtual VoxelType Get(int, int, int = 0, int = 0) const = 0;

  /// Get image value at an arbitrary discrete image location
  VoxelType operator ()(int) const;

  /// Get image value at an arbitrary discrete image location
  VoxelType operator ()(int, int, int = 0, int = 0) const;

  // Import overloaded non-virtual member functions from base class
  using ExtrapolateImageFunction::GetAsDouble;
  using ExtrapolateImageFunction::GetAsVector;

  /// Get scalar image value at an arbitrary discrete image location
  double GetAsDouble(int, int, int = 0, int = 0) const;

  /// Get vector image value at an arbitrary discrete image location
  void GetAsVector(Vector &, int, int, int = 0, int = 0) const;

};

////////////////////////////////////////////////////////////////////////////////
// Index extrapolation interface
////////////////////////////////////////////////////////////////////////////////

/**
 * Base class for generic voxel index extrapolation functions
 *
 * This abstract class is the base for extrapolate image functions which map
 * any given discrete voxel coordinate to one, which is inside the finite image
 * domain. The extrapolation function assigns the image value at the resulting
 * voxel coordinate to the initially mapped discrete voxel.
 *
 * This base class is templated over the type of the input image to be extrapolated
 * and thus subclasses can make use of image voxel type specific (non-virtual)
 * getters to access the image data. No conversion of the voxel type to a uniform
 * voxel type such as double scalar or vector type is required when the type of
 * the image to be extrapolated is known. Otherwise, if the image type is unknown,
 * use the abstract ExtrapolateImageFunction interface instead.
 *
 * \sa ExtrapolateImageFunction
 */
template <class TImage>
class IndexExtrapolateImageFunction
: public GenericExtrapolateImageFunction<TImage>
{
  mirtkAbstractMacro(IndexExtrapolateImageFunction);

  // ---------------------------------------------------------------------------
  // Types

public:

  typedef TImage                        ImageType; ///< Input image type
  typedef typename ImageType::VoxelType VoxelType; ///< Input voxel type
  typedef typename ImageType::RealType  RealType;  ///< Compatible floating-point type

  // ---------------------------------------------------------------------------
  // Attributes

protected:

  int _xmax; ///< Maximum index in x dimension
  int _ymax; ///< Maximum index in y dimension
  int _zmax; ///< Maximum index in z dimension
  int _tmax; ///< Maximum index in t dimension

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  IndexExtrapolateImageFunction();

public:

  /// Destructor
  virtual ~IndexExtrapolateImageFunction();

  // ---------------------------------------------------------------------------
  // Initialization

  /// Initialize image function
  virtual void Initialize();

  // ---------------------------------------------------------------------------
  // Voxel index transformation

  /// Transform voxel index such that it is inside the range [0, max]
  ///
  /// \param[in,out] idx Voxel index.
  /// \param[in]     max Maximum voxel index.
  virtual void TransformIndex(int &idx, int max) const = 0;

  /// Transform voxel index in x dimension
  void TransformX(int &) const;

  /// Transform voxel index in y dimension
  void TransformY(int &) const;

  /// Transform voxel index in z dimension
  void TransformZ(int &) const;

  /// Transform voxel index in t dimension
  virtual void TransformT(int &l) const;

  /// Transform voxel index in x dimension
  void TransformX(int &, int) const;

  /// Transform voxel index in y dimension
  void TransformY(int &, int) const;

  /// Transform voxel index in z dimension
  void TransformZ(int &, int) const;

  /// Transform voxel index in t dimension
  virtual void TransformT(int &l, int c) const;

  /// Transform 2D voxel index
  void Transform(int &, int &) const;

  /// Transform 3D voxel index
  void Transform(int &, int &, int &) const;

  /// Transform 4D voxel index
  void Transform(int &, int &, int &, int &) const;

  // ---------------------------------------------------------------------------
  // Extrapolation of discrete image function

  /// Whether voxel whose value is used is inside the foreground
  virtual bool IsForeground(int, int, int = 0, int = 0) const;

  // Import overloaded non-virtual member functions from base class
  using GenericExtrapolateImageFunction<TImage>::Get;

  /// Get image value at an arbitrary discrete image location
  virtual VoxelType Get(int, int, int = 0, int = 0) const;
  
};

////////////////////////////////////////////////////////////////////////////////
// Auxiliary macro for subclass implementation
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
#define mirtkExtrapolatorMacro(clsname, mode)                                  \
    mirtkObjectMacro(clsname);                                                 \
  public:                                                                      \
    /** Get extrapolation mode implemented by this extrapolator */             \
    inline virtual enum ExtrapolationMode ExtrapolationMode() const            \
    { return mode; }                                                           \
    /** Get extrapolation mode implemented by this class */                    \
    inline static  enum ExtrapolationMode ExtrapolationType()                  \
    { return mode; }                                                           \
  private:

// -----------------------------------------------------------------------------
#define mirtkGenericExtrapolatorMacro(clsname, mode)                           \
  mirtkExtrapolatorMacro(clsname, mode);                                       \
  public:                                                                      \
    /** Input image type */                                                    \
    typedef TImage                          ImageType;                         \
    /** Input voxel type */                                                    \
    typedef typename ImageType::VoxelType   VoxelType;                         \
    /** Compatible floating-point type */                                      \
    typedef typename ImageType::RealType    RealType;                          \
    /* Import overloaded non-virtual member functions from base class */       \
    using GenericExtrapolateImageFunction<TImage>::Get;                        \
  private:

////////////////////////////////////////////////////////////////////////////////
// Inline definitions -- ExtrapolateImageFunction
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline const ImageAttributes &ExtrapolateImageFunction::Attributes() const
{
  return Input()->Attributes();
}

// -----------------------------------------------------------------------------
inline int ExtrapolateImageFunction::X() const
{
  return Input()->X();
}

// -----------------------------------------------------------------------------
inline int ExtrapolateImageFunction::Y() const
{
  return Input()->Y();
}

// -----------------------------------------------------------------------------
inline int ExtrapolateImageFunction::Z() const
{
  return Input()->Z();
}

// -----------------------------------------------------------------------------
inline int ExtrapolateImageFunction::T() const
{
  return Input()->T();
}

// -----------------------------------------------------------------------------
inline int ExtrapolateImageFunction::N() const
{
  return Input()->N();
}

// -----------------------------------------------------------------------------
inline double ExtrapolateImageFunction::XSize() const
{
  return Input()->XSize();
}

// -----------------------------------------------------------------------------
inline double ExtrapolateImageFunction::YSize() const
{
  return Input()->YSize();
}

// -----------------------------------------------------------------------------
inline double ExtrapolateImageFunction::ZSize() const
{
  return Input()->ZSize();
}

// -----------------------------------------------------------------------------
inline double ExtrapolateImageFunction::TSize() const
{
  return Input()->TSize();
}

// -----------------------------------------------------------------------------
inline double ExtrapolateImageFunction
::Evaluate(double x, double y, double z, double t) const
{
  return this->GetAsDouble(static_cast<int>(round(x)),
                           static_cast<int>(round(y)),
                           static_cast<int>(round(z)),
                           static_cast<int>(round(t)));
}

// -----------------------------------------------------------------------------
inline double ExtrapolateImageFunction
::Evaluate(double x, double y, double z, double t)
{
  return const_cast<const ExtrapolateImageFunction *>(this)->Evaluate(x, y, z, t);
}

// -----------------------------------------------------------------------------
inline double ExtrapolateImageFunction::GetAsDouble(int i, int j, int k, int l) const
{
  // Must be implemented by subclass
  return this->_DefaultValue;
}

// -----------------------------------------------------------------------------
inline double ExtrapolateImageFunction::GetAsDouble(int idx) const
{
  int i, j, k, l;
  Input()->IndexToVoxel(idx, i, j, k, l);
  return this->GetAsDouble(i, j, k, l);
}

// -----------------------------------------------------------------------------
inline void ExtrapolateImageFunction
::GetAsVector(Vector &v, int i, int j, int k, int l) const
{
  // Must be implemented by subclass
  v = Vector(N(), this->_DefaultValue);
}

// -----------------------------------------------------------------------------
inline Vector ExtrapolateImageFunction::GetAsVector(int i, int j, int k, int l) const
{
  Vector v;
  this->GetAsVector(v, i, j, k, l);
  return v;
}

// -----------------------------------------------------------------------------
inline void ExtrapolateImageFunction::GetAsVector(Vector &v, int idx) const
{
  int i, j, k, l;
  Input()->IndexToVoxel(idx, i, j, k, l);
  this->GetAsVector(v, i, j, k, l);
}

// -----------------------------------------------------------------------------
inline Vector ExtrapolateImageFunction::GetAsVector(int idx) const
{
  Vector v;
  this->GetAsVector(v, idx);
  return v;
}

// -----------------------------------------------------------------------------
inline bool ExtrapolateImageFunction::IsForeground(int i, int j, int k, int l) const
{
  if (!this->Input()->HasBackgroundValue()) return true;
  const double value = this->GetAsDouble(i, j, k, l);
  const double bg    = this->Input()->GetBackgroundValueAsDouble();
  return (value != bg && (!IsNaN(value) || !IsNaN(bg)));
}

////////////////////////////////////////////////////////////////////////////////
// Inline definitions -- GenericExtrapolateImageFunction
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <class TImage>
inline GenericExtrapolateImageFunction<TImage>
::GenericExtrapolateImageFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
inline GenericExtrapolateImageFunction<TImage>
::~GenericExtrapolateImageFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void GenericExtrapolateImageFunction<TImage>
::Input(const BaseImage *input)
{
  ExtrapolateImageFunction::Input(dynamic_cast<const ImageType *>(input));
  if (input && !this->_Input) {
    cerr << this->NameOfClass() << "::Input: Invalid input image type" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline const TImage *GenericExtrapolateImageFunction<TImage>::Input() const
{
  return reinterpret_cast<const TImage *>(this->_Input);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void GenericExtrapolateImageFunction<TImage>::Initialize()
{
  // Initialize base class
  ExtrapolateImageFunction::Initialize();
  // Check type of input image -- also done by Input(const BaseImage *),
  //                              but just in case SetInput() has been used.
  if (!dynamic_cast<const ImageType *>(this->_Input)) {
    cerr << this->NameOfClass() << "::Initialize: Invalid input image type" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericExtrapolateImageFunction<TImage>::VoxelType
GenericExtrapolateImageFunction<TImage>::Get(int i, int j, int k, int l) const
{
  // Must be implemented by subclass
  return voxel_cast<VoxelType>(this->_DefaultValue);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericExtrapolateImageFunction<TImage>::VoxelType
GenericExtrapolateImageFunction<TImage>::Get(int idx) const
{
  int i, j, k, l;
  Input()->IndexToVoxel(idx, i, j, k, l);
  return this->Get(i, j, k, l);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericExtrapolateImageFunction<TImage>::VoxelType
GenericExtrapolateImageFunction<TImage>::operator ()(int idx) const
{
  return this->Get(idx);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericExtrapolateImageFunction<TImage>::VoxelType
GenericExtrapolateImageFunction<TImage>
::operator ()(int i, int j, int k, int l) const
{
  return this->Get(i, j, k, l);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline double GenericExtrapolateImageFunction<TImage>
::GetAsDouble(int i, int j, int k, int l) const
{
  return voxel_cast<double>(this->Get(i, j, k, l));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void GenericExtrapolateImageFunction<TImage>
::GetAsVector(Vector &v, int i, int j, int k, int l) const
{
  v = voxel_cast<Vector>(this->Get(i, j, k, l));
}

////////////////////////////////////////////////////////////////////////////////
// Inline definitions -- IndexExtrapolateImageFunction
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <class TImage>
IndexExtrapolateImageFunction<TImage>::IndexExtrapolateImageFunction()
:
  _xmax(0), _ymax(0), _zmax(0), _tmax(0)
{
}

// -----------------------------------------------------------------------------
template <class TImage>
IndexExtrapolateImageFunction<TImage>::~IndexExtrapolateImageFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
void IndexExtrapolateImageFunction<TImage>::Initialize()
{
  // Initialize base class
  GenericExtrapolateImageFunction<TImage>::Initialize();
  // Get range of indices
  _xmax = this->X() - 1;
  _ymax = this->Y() - 1;
  _zmax = this->Z() - 1;
  _tmax = this->T() - 1;
}

// -----------------------------------------------------------------------------
template <class TImage>
void IndexExtrapolateImageFunction<TImage>::TransformX(int &i) const
{
  this->TransformIndex(i, _xmax);
}

// -----------------------------------------------------------------------------
template <class TImage>
void IndexExtrapolateImageFunction<TImage>::TransformY(int &j) const
{
  this->TransformIndex(j, _ymax);
}

// -----------------------------------------------------------------------------
template <class TImage>
void IndexExtrapolateImageFunction<TImage>::TransformZ(int &k) const
{
  this->TransformIndex(k, _zmax);
}

// -----------------------------------------------------------------------------
template <class TImage>
void IndexExtrapolateImageFunction<TImage>::TransformT(int &l) const
{
  this->TransformIndex(l, _tmax);
}

// -----------------------------------------------------------------------------
template <class TImage>
void IndexExtrapolateImageFunction<TImage>::TransformX(int &i, int c) const
{
  i = c;
  this->TransformIndex(i, _xmax);
}

// -----------------------------------------------------------------------------
template <class TImage>
void IndexExtrapolateImageFunction<TImage>::TransformY(int &j, int c) const
{
  j = c;
  this->TransformIndex(j, _ymax);
}

// -----------------------------------------------------------------------------
template <class TImage>
void IndexExtrapolateImageFunction<TImage>::TransformZ(int &k, int c) const
{
  k = c;
  this->TransformIndex(k, _zmax);
}

// -----------------------------------------------------------------------------
template <class TImage>
void IndexExtrapolateImageFunction<TImage>::TransformT(int &l, int c) const
{
  l = c;
  this->TransformIndex(l, _tmax);
}

// -----------------------------------------------------------------------------
template <class TImage>
void IndexExtrapolateImageFunction<TImage>::Transform(int &i, int &j) const
{
  TransformX(i);
  TransformY(j);
}

// -----------------------------------------------------------------------------
template <class TImage>
void IndexExtrapolateImageFunction<TImage>::Transform(int &i, int &j, int &k) const
{
  TransformX(i);
  TransformY(j);
  TransformZ(k);
}

// -----------------------------------------------------------------------------
template <class TImage>
void IndexExtrapolateImageFunction<TImage>::Transform(int &i, int &j, int &k, int &l) const
{
  TransformX(i);
  TransformY(j);
  TransformZ(k);
  TransformT(l);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline bool IndexExtrapolateImageFunction<TImage>
::IsForeground(int i, int j, int k, int l) const
{
  Transform(i, j, k, l);
  return this->Input()->IsForeground(i, j, k, l);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename IndexExtrapolateImageFunction<TImage>::VoxelType
IndexExtrapolateImageFunction<TImage>::Get(int i, int j, int k, int l) const
{
  Transform(i, j, k, l);
  return this->Input()->Get(i, j, k, l);
}


} // namespace mirtk

////////////////////////////////////////////////////////////////////////////////
// Extrapolation functions
////////////////////////////////////////////////////////////////////////////////

#include "mirtk/ConstExtrapolateImageFunction.h"
#include "mirtk/ConstExtrapolateImageFunctionWithPeriodicTime.h"
#include "mirtk/NearestNeighborExtrapolateImageFunction.h"
#include "mirtk/RepeatExtrapolateImageFunction.h"
#include "mirtk/MirrorExtrapolateImageFunction.h"

////////////////////////////////////////////////////////////////////////////////
// Instantiation
////////////////////////////////////////////////////////////////////////////////

namespace mirtk {


// -----------------------------------------------------------------------------
template <class TImage>
GenericExtrapolateImageFunction<TImage> *
GenericExtrapolateImageFunction<TImage>
::New(enum ExtrapolationMode mode, const TImage *image)
{
  GenericExtrapolateImageFunction<TImage> *p = NULL;
  switch (mode) {
    case Extrapolation_None:   { p = NULL;                                                         break; }
    case Extrapolation_Const:  { p = new GenericConstExtrapolateImageFunction          <TImage>(); break; }
    case Extrapolation_NN:     { p = new GenericNearestNeighborExtrapolateImageFunction<TImage>(); break; }
    case Extrapolation_Repeat: { p = new GenericRepeatExtrapolateImageFunction         <TImage>(); break; }
    case Extrapolation_Mirror: { p = new GenericMirrorExtrapolateImageFunction         <TImage>(); break; }
    case Extrapolation_ConstWithPeriodicTime:
      p = new GenericConstExtrapolateImageFunctionWithPeriodicTime<TImage>();
      break;
    default:
   	  cerr << "GenericExtrapolateImageFunction::New: Unknown extrapolation mode: " << mode << endl;
      exit(1);
  }
  if (p) p->Input(image);
  return p;
}


} // namespace mirtk

#endif // MIRTK_ExtrapolateImageFunction_H
