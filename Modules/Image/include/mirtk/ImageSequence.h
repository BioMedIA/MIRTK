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

#ifndef MIRTK_ImageSequence_H
#define MIRTK_ImageSequence_H

#include "mirtk/Object.h"
#include "mirtk/BaseImage.h"
#include "mirtk/Array.h"


namespace mirtk {


////////////////////////////////////////////////////////////////////////////////
// Image channel
////////////////////////////////////////////////////////////////////////////////

/**
 * Single channel of a multi-channel image
 *
 * \note The image associated with the channel cannot have a fourth dimension!
 */
template <class TImage = BaseImage>
class ImageChannel : public Object
{
  mirtkObjectMacro(ImageChannel);

public:

  // ---------------------------------------------------------------------------
  // Types

  /// Type of image associated with this channel
  typedef TImage   ImageType;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Constructor
  ///
  /// \param[in] image  Image associated with channel.
  /// \param[in] manage Whether to manage memory of image.
  /// \param[in] copy   Whether to copy image if \p manage is \c false.
  ImageChannel(ImageType *image = NULL, bool manage = false, bool copy = false);

  /// Copy constructor
  ///
  /// \attention If the given channel object owns the associated image and thus
  ///            would destruct it upon destruction of the channel itself, this
  ///            copy constructor takes over ownership of the channel image.
  ///            This is required for using ImageChannel as type for an STL
  ///            container which creates copies of the elements.
  ImageChannel(const ImageChannel &);

  /// Copy constructor
  ///
  /// \attention If the given channel object owns the associated image and thus
  ///            would destruct it upon destruction of the channel itself, this
  ///            copy constructor takes over ownership of the channel image.
  ///            This is required for using ImageChannel as type for an STL
  ///            container which creates copies of the elements.
  ///
  /// \note Implemented only for TImage template argument BaseImage.
  template <class TOtherImage>
  ImageChannel(const ImageChannel<TOtherImage> &);

  /// Assignment operator
  ///
  /// \attention If the given channel object owns the associated image and thus
  ///            would destruct it upon destruction of the channel itself, this
  ///            assignment transfer also the ownership of the channel image.
  ///            This is required for using ImageChannel as type for an STL
  ///            container which creates copies of the elements.
  ImageChannel &operator =(const ImageChannel &);

  /// Assignment operator
  ///
  /// \attention If the given channel object owns the associated image and thus
  ///            would destruct it upon destruction of the channel itself, this
  ///            assignment transfer also the ownership of the channel image.
  ///            This is required for using ImageChannel as type for an STL
  ///            container which creates copies of the elements.
  ///
  /// \note Implemented only for TImage template argument BaseImage.
  template <class TOtherImage>
  ImageChannel &operator =(const ImageChannel<TOtherImage> &);

  /// Destructor
  virtual ~ImageChannel();

  // ---------------------------------------------------------------------------
  // Channel image

  /// Set image associated with this channel
  ///
  /// \param[in] image  Image associated with this channel.
  /// \param[in] manage Whether to manage memory of image.
  /// \param[in] copy   Whether to copy image if \p manage is \c false.
  void Image(ImageType *image, bool manage = false, bool copy = false);

  /// Get image associated with this channel
  ImageType *Image() const;

  // ---------------------------------------------------------------------------
  // Data members
protected:

  ImageType    *_Image;  ///< Image associated with this channel
  mutable bool  _Manage; ///< Whether to manage the associated image

  // Copy constructor and assignment operator modify the _ManageMemory
  // member of the const channel object which is being copied
  template <class T> friend class ImageChannel;
};

////////////////////////////////////////////////////////////////////////////////
// Multi-channel image
////////////////////////////////////////////////////////////////////////////////

/**
 * Auxiliary type to store channels of a single frame of a temporal image sequence
 *
 * All channels of a given frame have the same spatial image attributes. Moreover,
 * the temporal origin and voxel size must be identical for all channels as well.
 * The size of the channel images in the time domain is always one.
 *
 * \note This data type can also be used for single multi-channel images which are
 *       not part of a temporal sequence of images.
 */
template <class TChannel = ImageChannel<BaseImage> >
class ImageFrame : public Object
{
  mirtkObjectMacro(ImageFrame);

public:

  // ---------------------------------------------------------------------------
  // Types

  typedef typename TChannel::ImageType ImageType;   ///< Type of channel image
  typedef TChannel                     ChannelType; ///< Type of a channel

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  ImageFrame();

  /// Copy constructor
  ImageFrame(const ImageFrame &);

  /// Assignment operator
  ImageFrame &operator =(const ImageFrame &);

  /// Destructor
  virtual ~ImageFrame();

  // ---------------------------------------------------------------------------
  // Frame

  /// Add image as channel to this frame
  ///
  /// \param[in] image  Channel image. If the image has a fourth dimension, a
  ///                   copy of each 2D/3D sub-image is added as separate channel.
  /// \param[in] copy   Whether to copy image even if it is 2D/3D only. Otherwise,
  ///                   the given image is associated with the channel which also
  ///                   takes over its ownership, i.e., destructing it when the
  ///                   channel object itself is destructed.
  /// \param[in] manage Whether to manage memory of input image if not copied.
  void Add(ImageType *image, bool manage = false, bool copy = false);

  /// Clear frame
  void Clear();

  // ---------------------------------------------------------------------------
  // Channels

  /// Set number of image channels
  ///
  /// Newly added channels are initialized using the default constructor of the
  /// ChannelType, i.e., the associated image pointers set to NULL. Existing
  /// channels within the specified range of channels are not changed.
  void NumberOfChannels(int);

  /// Get number of image channels
  int NumberOfChannels() const;

  /// Get image channel
  ChannelType &Channel(int = 0);

  /// Get image channel
  const ChannelType &Channel(int = 0) const;

  /// Set channel image
  ///
  /// \param[in] idx    Index of channel.
  /// \param[in] image  Image to associate with specified channel.
  /// \param[in] manage Whether to manage memory of input image.
  /// \param[in] copy   Whether to copy image if \c manage is \c false.
  void Image(int idx, ImageType *image, bool manage = false, bool copy = false);

  /// Get channel image
  ImageType *Image(int = 0) const;

  // ---------------------------------------------------------------------------
  // Attributes

  /// Get attributes of this frame
  ///
  /// The spatial attributes of the frame correspond to those of its channels
  /// which are the same for all channels.
  ///
  /// The temporal attributes are:
  /// - _t:       Number of channels
  /// - _dt:      Temporal voxel size of frame (in ms)
  /// - _torigin: Time coordinate of frame (in ms)
  ImageAttributes Attributes() const;

  /// Get number of voxels per channel
  int NumberOfVoxels() const;

  /// Get number of voxels in x dimension for each channel
  int X() const;

  /// Get number of voxels in y dimension for each channel
  int Y() const;

  /// Get number of voxels in z dimension for each channel
  int Z() const;

  /// Get size of each voxel in x dimension
  double XSize() const;

  /// Get size of each voxel in y dimension
  double YSize() const;

  /// Get size of each voxel in z dimension
  double ZSize() const;

  /// Get voxel size
  void GetPixelSize(double *, double *, double *) const;

  /// Convert voxel to world coordinate
  void ImageToWorld(double &, double &, double &) const;

  /// Convert world to voxel coordinate
  void WorldToImage(double &, double &, double &) const;

  /// Get time of temporal frame
  double Time() const;

  // ---------------------------------------------------------------------------
  // Attributes
private:

  Array<ChannelType> _Channel; ///< Channels of image frame

};

////////////////////////////////////////////////////////////////////////////////
// Temporal multi-channel image sequence
////////////////////////////////////////////////////////////////////////////////

/**
 * Auxiliary type to store ordered sequence of temporal (multi-channel) image frames
 */
template <class TFrame = ImageFrame<> >
class ImageSequence : public Object
{
  mirtkObjectMacro(ImageSequence);

public:

  // ---------------------------------------------------------------------------
  // Types

  typedef typename TFrame::ImageType     ImageType;   ///< Type of images
  typedef typename TFrame::ChannelType   ChannelType; ///< Type of channels
  typedef TFrame                         FrameType;   ///< Type of frames

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Constructor
  ImageSequence();

  /// Copy constructor
  ImageSequence(const ImageSequence &);

  /// Assignment operator
  ImageSequence &operator =(const ImageSequence &);

  /// Destructor
  ~ImageSequence();

  // ---------------------------------------------------------------------------
  // Sequence

  /// Add channels or frames to image sequence
  ///
  /// \param[in] image  Image to add to this sequence. If the image has a fouth
  ///                   dimension, a copy of each 2D/3D sub-image is added.
  /// \param[in] manage Whether to manage memory of input image if not copied.
  /// \param[in] copy   Whether to copy image if \p manage is \c false even
  ///                   though it has no fourth dimension.
  void Add(ImageType *image, bool manage = false, bool copy = false);

  /// Clear image sequence
  void Clear();

  // ---------------------------------------------------------------------------
  // Frames

  /// Set number of frames
  void NumberOfFrames(int) const;

  /// Get number of frames
  int NumberOfFrames() const;

  /// Get frame of sequence
  FrameType &Frame(int);

  /// Get frame of sequence
  const FrameType &Frame(int) const;

  // ---------------------------------------------------------------------------
  // Channels

  /// Set number of channels per frame
  void NumberOfChannels(int) const;

  /// Get number of channels per frame
  int NumberOfChannels() const;

  /// Get channel
  ///
  /// \param[in] idx Index of channel across all frames.
  ChannelType &Channel(int idx);

  /// Get channel
  ///
  /// \param[in] idx Index of channel across all frames.
  const ChannelType &Channel(int idx) const;

  /// Get channel
  ///
  /// \param[in] f Index of frame.
  /// \param[in] c Index of channel.
  ChannelType &Channel(int f, int c);

  /// Get channel
  ///
  /// \param[in] f Index of frame.
  /// \param[in] c Index of channel.
  const ChannelType &Channel(int f, int c) const;

  // ---------------------------------------------------------------------------
  // Images

  /// Get number of frames times number of channels
  int NumberOfImages() const;

  /// Set image associated with channel
  ///
  /// \param[in] f      Index of frame.
  /// \param[in] c      Index of channel.
  /// \param[in] image  Image which is associated with the channel.
  /// \param[in] manage Whether to manage the memory of the image.
  /// \param[in] copy   Whether to copy the image if \p manage is \c false.
  void Image(int f, int c, ImageType *image, bool manage, bool copy);

  /// Set image associated with channel to copy of given image
  ///
  /// \param[in] f     Index of frame.
  /// \param[in] c     Index of channel.
  /// \param[in] image Image of which a copy is associated with the channel.
  void Image(int f, int c, const ImageType *image);

  /// Get image associated with channel
  ///
  /// \param[in] idx Index of channel across all frames.
  ImageType *Image(int idx) const;

  /// Get image associated with channel
  ///
  /// \param[in] f Index of frame.
  /// \param[in] c Index of channel.
  ImageType *Image(int f, int c) const;

  // ---------------------------------------------------------------------------
  // Attributes

  /// Get attributes of this image sequence
  ///
  /// The spatial attributes of the sequence correspond to those of its frames
  /// which are the same for all frames. The temporal attributes are:
  /// - _t:       Number of frames
  /// - _dt:      Temporal size   of first frame (in ms)
  /// - _torigin: Time coordinate of first frame (in ms)
  ImageAttributes Attributes() const;

  /// Get number of voxels per channel
  int NumberOfVoxels() const;

  /// Get number of voxels in x dimension for each channel
  int X() const;

  /// Get number of voxels in y dimension for each channel
  int Y() const;

  /// Get number of voxels in z dimension for each channel
  int Z() const;

  /// Get number of frames
  int T() const;

  /// Get size of each voxel in x dimension
  double XSize() const;

  /// Get size of each voxel in y dimension
  double YSize() const;

  /// Get size of each voxel in z dimension
  double ZSize() const;

  /// Get voxel size
  void GetPixelSize(double *, double *, double *) const;

  /// Convert voxel to world coordinate
  void ImageToWorld(double &, double &, double &) const;

  /// Convert world to voxel coordinate
  void WorldToImage(double &, double &, double &) const;

  /// Get time of temporal frame
  double Time(int) const;

  // ---------------------------------------------------------------------------
  // Attributes
private:

  Array<FrameType> _Frame; ///< Frames of image frame

};


} // namespace mirtk


// Inline definitions
#include "mirtk/ImageSequence.hh"


#endif // MIRTK_ImageSequence_H
