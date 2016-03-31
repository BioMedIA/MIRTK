/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
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

#ifndef MIRTK_ImageToImage_H
#define MIRTK_ImageToImage_H

#include "mirtk/Object.h"
#include "mirtk/GenericImage.h"


namespace mirtk {


/**
 * Abstract base class for any general image to image filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take an image as input and produce an image as output. Each
 * derived class has to implement all abstract member functions.
 */
template <class TVoxel>
class ImageToImage : public Object
{
  mirtkAbstractMacro(ImageToImage);

  // ---------------------------------------------------------------------------
  // Types

public:

  /// Input/output image voxel type
  typedef TVoxel VoxelType;

  /// Input/output image type
  typedef GenericImage<VoxelType> ImageType;

  // ---------------------------------------------------------------------------
  // Attributes

protected:

  /// Input image for filter
  mirtkPublicAggregateMacro(const ImageType, Input);

  /// Output image for filter
  mirtkPublicAggregateMacro(ImageType, Output);

  /// Buffer
  mirtkAggregateMacro(ImageType, Buffer);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  ImageToImage();

  /// Destructor
  virtual ~ImageToImage();

  // ---------------------------------------------------------------------------
  // Evaluation

protected:

  /// Initialize the filter. This function must be called by any derived
  /// filter class to perform some initialize tasks.
  void Initialize(bool);

  /// Initialize the filter. This function must be called by any derived
  /// filter class to perform some initialize tasks.
  virtual void Initialize();

  /// Finalize the filter. This function must be called by any derived
  /// filter class to perform some initialize tasks.
  virtual void Finalize();

public:

  /// Run filter on entire image
  virtual void Run();

  /// Run filter on single voxel
  virtual double Run(int, int, int, int = 0);

  /// Returns whether the filter requires buffering. Any derived class must
  /// implement this member function to indicate whether the filter should
  /// buffer the input in case that input and output are equal. For example,
  /// filters which only require the voxel value to calculate their output
  /// should return false, otherwise true.
  virtual bool RequiresBuffering() const;

};

////////////////////////////////////////////////////////////////////////////////
// Auxiliary macro for subclass implementation
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
#define mirtkDefineImageFilterTypesMacro(voxeltype)                            \
  public:                                                                      \
    /** Type of image filter base class */                                     \
    typedef mirtk::ImageToImage<voxeltype> Baseclass;                          \
    /** Input/output image voxel type */                                       \
    typedef typename Baseclass::VoxelType VoxelType;                           \
    /** Input/output image type */                                             \
    typedef typename Baseclass::ImageType ImageType;                           \
  private:

// -----------------------------------------------------------------------------
#define mirtkAbstractImageFilterMacro(clsname, voxeltype)                      \
  mirtkAbstractMacro(clsname);                                                 \
  mirtkDefineImageFilterTypesMacro(voxeltype);                                 \
  public:                                                                      \
    typedef clsname<voxeltype> Superclass;                                     \
  private:

// -----------------------------------------------------------------------------
#define mirtkImageFilterMacro(clsname, voxeltype)                              \
  mirtkObjectMacro(clsname);                                                   \
  mirtkDefineImageFilterTypesMacro(voxeltype);                                 \
  protected:                                                                   \
    /** Indicates that this filter cannot run in-place */                      \
    virtual bool RequiresBuffering() const { return true; }                    \
  private:

// -----------------------------------------------------------------------------
#define mirtkInPlaceImageFilterMacro(clsname, voxeltype)                       \
  mirtkObjectMacro(clsname);                                                   \
  mirtkDefineImageFilterTypesMacro(voxeltype);                                 \
  protected:                                                                   \
    /** Indicates that this filter can run in-place */                         \
    virtual bool RequiresBuffering() const { return false; }                   \
  private:


} // namespace mirtk

#endif // MIRTK_ImageToImage_H
