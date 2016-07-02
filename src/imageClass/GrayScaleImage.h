#pragma once

#include <string>
#include <vector>

#include <cassert>

//*******************************************************************************************************************
/*!
* A gray scale image
*
* Data can be loaded and stored from/to png files
*/
//*******************************************************************************************************************
class GrayScaleImage
{
public:
  typedef unsigned char pixel_t;

  /// Creates an empty image: width and height are the amount of pixels in x/y coordinate
  GrayScaleImage( unsigned int _width, unsigned int _height );
  
  /// Load image from a grayscale png file
  GrayScaleImage( const std::string & pngFilename );

  /// Save current image to png file
  void save( const std::string & pngFilename );

  /// Returns a resized version the image.
  GrayScaleImage getResizedImage( unsigned int newWidth, unsigned int newHeight, bool bilinear=true ) const;

  /// Number of pixels in x direction
  unsigned int width()  const { return size_[0]; }
  /// Number of pixels in y direction
  unsigned int height() const { return size_[1]; }

  /// size(0) == width(), size(1) == height()
  unsigned int size( unsigned int coord ) const;

  /// Returns the pixel value at position x,y. Returned values are between 0.0 (white) and 1.0 (black)
  double operator() ( int x, int y ) const;

  /// Sets the value at the given pixel. 'val' has to be between 0 and 1
  void setElement( int x, int y, double val);

  /// Access raw pixel values - values are between 0 (black) and 255(white)
  unsigned char & getElement ( int x, int y );
  unsigned char   getElement ( int x, int y ) const;


  static pixel_t pixelValueFromString( const std::string & str );

protected:
  GrayScaleImage() {}

  unsigned int size_[2];                   //< 0=width,  1=height
  std::vector<unsigned char> image_; //< raw pixels
};



inline unsigned char &  GrayScaleImage::getElement ( int x, int y ) {
  assert( x >= 0 && x < int(size_[0]) );
  assert( y >= 0 && y < int(size_[1]) );
  const unsigned int yFlip = size_[1] - unsigned(y) - unsigned(1);
  return image_[ yFlip * size_[0] + unsigned(x) ];
}

inline unsigned char GrayScaleImage::getElement ( int x, int y ) const {
  return const_cast<GrayScaleImage*> ( this )->getElement(x,y);
}


inline unsigned int GrayScaleImage::size( unsigned int coord ) const
{
  assert( coord < 2);
  return size_[coord];
}




