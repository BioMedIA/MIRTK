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

#ifndef MIRTK_ImageSequence_HH
#define MIRTK_ImageSequence_HH

#include "mirtk/ImageSequence.h"
#include "mirtk/Stream.h"


namespace mirtk {


////////////////////////////////////////////////////////////////////////////////
// Inline definitions of ImageChannel
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <>
ImageChannel<irtkBaseImage>::ImageChannel(ImageType *image, bool manage, bool copy)
:
  _Image (image),
  _Manage(manage)
{
  if (image) {
    if (image->GetT() > 0) {
      cerr << "ImageChannel::ImageChannel: Channel image cannot have fourth dimension" << endl;
      exit(1);
    }
    if (copy && !manage) {
      _Image  = BaseImage::New(image);
      _Manage = true;
    }
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
ImageChannel<TImage>::ImageChannel(ImageType *image, bool manage, bool copy)
:
  _Image (image),
  _Manage(manage)
{
  if (image) {
    if (image->GetT() > 0) {
      cerr << "ImageChannel::ImageChannel: Channel image cannot have fourth dimension" << endl;
      exit(1);
    }
    if (copy && !manage) {
      _Image  = new ImageType(image);
      _Manage = true;
    }
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
ImageChannel<TImage> &ImageChannel<TImage>::operator =(const ImageChannel<TImage> &other)
{
  if (_Manage) delete _Image;
  _Image  = other._Image;
  _Manage = other._Manage;
  other._Manage = false; // take over ownership
  return *this;
}

// -----------------------------------------------------------------------------
template <class TImage>
ImageChannel<TImage>::ImageChannel(const ImageChannel<TImage> &other)
:
  Object(other),
  _Image (NULL),
  _Manage(false)
{
  *this = other;
}

// -----------------------------------------------------------------------------
template <> template <class TOtherImage>
ImageChannel<irtkBaseImage> &ImageChannel<irtkBaseImage>::operator =(const ImageChannel<TOtherImage> &other)
{
  if (_Manage) delete _Image;
  _Image  = other._Image;
  _Manage = other._Manage;
  other._Manage = false; // take over ownership
  return *this;
}

// -----------------------------------------------------------------------------
template <> template <class TOtherImage>
ImageChannel<irtkBaseImage>::ImageChannel(const ImageChannel<TOtherImage> &other)
:
  Object(other),
  _Image (NULL),
  _Manage(false)
{
  *this = other;
}

// -----------------------------------------------------------------------------
template <class TImage>
ImageChannel<TImage>::~ImageChannel()
{
  if (_Manage) delete _Image;
}

// -----------------------------------------------------------------------------
template <>
void ImageChannel<irtkBaseImage>::Image(ImageType *image, bool manage, bool copy)
{
  if (_Manage) delete _Image;
  if (image && copy && !manage) {
    _Image  = BaseImage::New(image);
    _Manage = true;
  } else {
    _Image  = image;
    _Manage = manage;
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
void ImageChannel<TImage>::Image(ImageType *image, bool manage, bool copy)
{
  if (_Manage) delete _Image;
  if (image && copy && !manage) {
    _Image  = new ImageType(*image);
    _Manage = true;
  } else {
    _Image  = image;
    _Manage = manage;
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
TImage *ImageChannel<TImage>::Image() const
{
  return _Image;
}

////////////////////////////////////////////////////////////////////////////////
// Inline definitions of ImageFrame
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class TChannel>
inline ImageFrame<TChannel>::ImageFrame()
{
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline ImageFrame<TChannel>::~ImageFrame()
{
}

// =============================================================================
// Frame
// =============================================================================

// -----------------------------------------------------------------------------
template <class TChannel>
inline void ImageFrame<TChannel>::Add(ImageType *image, bool manage, bool copy)
{
  if (_Channel.size() > 0) {
    ImageAttributes attr = image->Attributes();
    attr._t  = 1; // may differ
    attr._dt =  _Channel[0].Image().GetTSize();
    if (attr != _Channel[0].Image().Attributes()) {
      cerr << "ImageFrame::Add: Attributes of image do not match those of first channel" << endl;
      exit(1);
    }
  }
  const int T = image->GetT();
  if (T > 0) {
    for (int t = 0; t < T; ++t) {
      BaseImage *channel;
      image.GetFrame(t, channel);
      _Channel.push_back(ChannelType(channel, true, false));
    }
  } else {
    _Channel.push_back(ChannelType(image, manage, copy));
  }
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline void ImageFrame<TChannel>::Clear()
{
  _Channel.clear();
}

// =============================================================================
// Channels
// =============================================================================

// -----------------------------------------------------------------------------
template <class TChannel>
inline void ImageFrame<TChannel>::NumberOfChannels(int n)
{
  _Channel.resize(n);
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline int ImageFrame<TChannel>::NumberOfChannels() const
{
  return static_cast<int>(_Channel.size());
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline typename ImageFrame<TChannel>::ChannelType &ImageFrame<TChannel>::Channel(int idx)
{
  return _Channel[idx];
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline const typename ImageFrame<TChannel>::ChannelType &ImageFrame<TChannel>::Channel(int idx) const
{
  return _Channel[idx];
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline void ImageFrame<TChannel>::Image(int idx, ImageType *image, bool manage, bool copy)
{
  _Channel[idx].Image(image, manage, copy);
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline typename ImageFrame<TChannel>::ImageType *ImageFrame<TChannel>::Image(int idx) const
{
  return _Channel[idx].Image();
}

// =============================================================================
// Attributes
// =============================================================================

// -----------------------------------------------------------------------------
template <class TChannel>
inline ImageAttributes ImageFrame<TChannel>::Attributes() const
{
  const ImageType *image = Image(0);
  return image ? image->GetImageAttributes() : ImageAttributes();
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline int ImageFrame<TChannel>::NumberOfVoxels() const
{
  const ImageType *image = Image(0);
  return image ? image->GetNumberOfVoxels() : 0;
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline int ImageFrame<TChannel>::X() const
{
  const ImageType *image = Image(0);
  return image ? image->GetX() : 0;
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline int ImageFrame<TChannel>::Y() const
{
  const ImageType *image = Image(0);
  return image ? image->GetY() : 0;
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline int ImageFrame<TChannel>::Z() const
{
  const ImageType *image = Image(0);
  return image ? image->GetZ() : 0;
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline double ImageFrame<TChannel>::XSize() const
{
  const ImageType *image = Image(0);
  return image ? image->GetXSize() : .0;
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline double ImageFrame<TChannel>::YSize() const
{
  const ImageType *image = Image(0);
  return image ? image->GetYSize() : .0;
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline double ImageFrame<TChannel>::ZSize() const
{
  const ImageType *image = Image(0);
  return image ? image->GetZSize() : .0;
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline void ImageFrame<TChannel>::GetPixelSize(double *dx, double *dy, double *dz) const
{
  const ImageType *image = Image(0);
  if (image) image->GetPixelSize(dx, dy, dz);
  else {
    if (dx) *dx = .0;
    if (dy) *dy = .0;
    if (dz) *dz = .0;
  }
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline void ImageFrame<TChannel>::ImageToWorld(double &x, double &y, double &z) const
{
  const ImageType *image = Image(0);
  if (image) image->ImageToWorld(x, y, z);
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline void ImageFrame<TChannel>::WorldToImage(double &x, double &y, double &z) const
{
  const ImageType *image = Image(0);
  if (image) image->WorldToImage(x, y, z);
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline double ImageFrame<TChannel>::Time() const
{
  const ImageType *image = Image(0);
  if (image) image->ImageToTime(.0);
}

////////////////////////////////////////////////////////////////////////////////
// Inline definitions of ImageSequence
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class TFrame>
inline ImageSequence<TFrame>::ImageSequence()
{
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline ImageSequence<TFrame>::ImageSequence(const ImageSequence &other)
:
  Object(other),
  _Frame(other._Frame)
{
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline ImageSequence<TFrame> &ImageSequence<TFrame>::operator =(const ImageSequence &other)
{
  _Frame = other._Frame;
  return *this;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline ImageSequence<TFrame>::~ImageSequence()
{
}

// =============================================================================
// Sequence
// =============================================================================

// -----------------------------------------------------------------------------
template <class TFrame>
inline void ImageSequence<TFrame>::Add(ImageType *image, bool manage, bool copy)
{
  if (NumberOfFrames() > 0 && NumberOfChannels() > 0) {
    if (!image->HasSpatialAttributesOf(Image(0, 0))) {
      cerr << "ImageSequence::Add: Spatial attributes of image differ from those of first frame" << endl;
      exit(1);
    }
  }
  // Reserve enough entries in vector to ensure that no reallocation
  // takes place during the insert to keep tmp iterator below valid
  _Frame.reserve(_Frame.size() + image->GetT());
  // Insert each input frame (dt != 0) or channel (dt == 0)
  // Frames are sorted by increasing time and channels appended
  for (int l = 0; l < image->GetT(); ++l) {
    // Find corresponding frame or add new one if necessary
    const double                              time  = image->ImageToTime(l);
    typename Array<FrameType>::iterator frame = _Frame.begin();
    while (frame != _Frame.end()) {
      const double t = frame->Time();
      if      (t == time) break;
      else if (t  > time) {
        typename Array<FrameType>::iterator tmp = frame - 1;
        _Frame.insert(frame, FrameType());
        frame = ++tmp;
        break;
      }
      ++frame;
    }
    if (frame == _Frame.end()) {
      _Frame.push_back(FrameType());
      frame = _Frame.end() - 1;
    }
    // Add channel to frame
    if (image->GetT() > 0) {
      BaseImage *channel;
      image->GetFrame(l, channel);
      frame->Add(channel, true, false);
    } else {
      frame->Add(image, manage, copy);
    }
  }
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline void ImageSequence<TFrame>::Clear()
{
  _Frame.clear();
}

// =============================================================================
// Frames
// =============================================================================

// -----------------------------------------------------------------------------
template <class TFrame>
inline void ImageSequence<TFrame>::NumberOfFrames(int n) const
{
  _Frame.resize(n);
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline int ImageSequence<TFrame>::NumberOfFrames() const
{
  return static_cast<int>(_Frame.size());
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline typename ImageSequence<TFrame>::FrameType &ImageSequence<TFrame>::Frame(int f)
{
  return _Frame[f];
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline const typename ImageSequence<TFrame>::FrameType &ImageSequence<TFrame>::Frame(int f) const
{
  return _Frame[f];
}

// =============================================================================
// Channels
// =============================================================================

// -----------------------------------------------------------------------------
template <class TFrame>
inline void ImageSequence<TFrame>::NumberOfChannels(int n) const
{
  for (int f = 0; f < NumberOfChannels(); ++f) Frame(f).NumberOfChannels(n);
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline int ImageSequence<TFrame>::NumberOfChannels() const
{
  return NumberOfFrames() > 0 ? Frame(0).NumberOfChannels() : 0;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline typename ImageSequence<TFrame>::ChannelType &ImageSequence<TFrame>::Channel(int idx)
{
  const int num = NumberOfChannels();
  return _Frame[idx / num][idx % num];
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline const typename ImageSequence<TFrame>::ChannelType &ImageSequence<TFrame>::Channel(int idx) const
{
  const int num = NumberOfChannels();
  return _Frame[idx / num][idx % num];
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline typename ImageSequence<TFrame>::ChannelType &ImageSequence<TFrame>::Channel(int f, int c)
{
  return _Frame[f][c];
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline const typename ImageSequence<TFrame>::ChannelType &ImageSequence<TFrame>::Channel(int f, int c) const
{
  return _Frame[f][c];
}

// =============================================================================
// Images
// =============================================================================

// -----------------------------------------------------------------------------
template <class TFrame>
inline int ImageSequence<TFrame>::NumberOfImages() const
{
  return NumberOfFrames() * NumberOfChannels();
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline typename ImageSequence<TFrame>::ImageType *ImageSequence<TFrame>::Image(int idx) const
{
  return Channel(idx).Image();
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline typename ImageSequence<TFrame>::ImageType *ImageSequence<TFrame>::Image(int f, int c) const
{
  return Channel(f, c).Image();
}

// =============================================================================
// Attributes
// =============================================================================

// -----------------------------------------------------------------------------
template <class TFrame>
inline ImageAttributes ImageSequence<TFrame>::Attributes() const
{
  const ImageType *image = Image(0, 0);
  ImageAttributes attr;
  if (image) attr = image->GetImageAttributes();
  attr._t = NumberOfFrames();
  return attr;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline int ImageSequence<TFrame>::NumberOfVoxels() const
{
  const ImageType *image = Image(0, 0);
  return image ? image->NumberOfVoxels() : 0;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline int ImageSequence<TFrame>::X() const
{
  const ImageType *image = Image(0, 0);
  return image ? image->GetX() : 0;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline int ImageSequence<TFrame>::Y() const
{
  const ImageType *image = Image(0, 0);
  return image ? image->GetY() : 0;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline int ImageSequence<TFrame>::Z() const
{
  const ImageType *image = Image(0, 0);
  return image ? image->GetZ() : 0;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline int ImageSequence<TFrame>::T() const
{
  return NumberOfFrames();
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline double ImageSequence<TFrame>::XSize() const
{
  const ImageType *image = Image(0, 0);
  return image ? image->GetXSize() : .0;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline double ImageSequence<TFrame>::YSize() const
{
  const ImageType *image = Image(0, 0);
  return image ? image->GetYSize() : .0;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline double ImageSequence<TFrame>::ZSize() const
{
  const ImageType *image = Image(0, 0);
  return image ? image->GetZSize() : .0;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline void ImageSequence<TFrame>::GetPixelSize(double *dx, double *dy, double *dz) const
{
  const ImageType *image = Image(0, 0);
  if (image) image->GetPixelSize(dx, dy, dz);
  else {
    if (dx) *dx = .0;
    if (dy) *dy = .0;
    if (dz) *dz = .0;
  }
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline void ImageSequence<TFrame>::ImageToWorld(double &x, double &y, double &z) const
{
  const ImageType *image = Image(0, 0);
  if (image) image->ImageToWorld(x, y, z);
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline void ImageSequence<TFrame>::WorldToImage(double &x, double &y, double &z) const
{
  const ImageType *image = Image(0, 0);
  if (image) image->WorldToImage(x, y, z);
}


} // namespace mirtk

#endif // MIRTK_ImageSequence_HH
