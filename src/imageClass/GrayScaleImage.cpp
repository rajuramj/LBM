#include "GrayScaleImage.h"
#include "lodepng.h"

#include <sstream>
#include <stdexcept>
#include <limits>

GrayScaleImage::GrayScaleImage( unsigned int _width, unsigned int _height )
{
    size_[0] = _width;
    size_[1] = _height;
    image_.resize( size_[0] * size_[1] );
}

GrayScaleImage::GrayScaleImage( const std::string & pngFilename )
{
    unsigned int tmpWidth, tmpHeight;
    unsigned int error = lodepng::decode( image_, tmpWidth, tmpHeight, pngFilename, LCT_GREY, 8 );
    size_[0] = tmpWidth;
    size_[1] = tmpHeight;

    if ( error )
      throw std::invalid_argument( std::string( "Error while loading PNG file: " ) + lodepng_error_text(error)  );
}

void GrayScaleImage::save( const std::string & pngFilename )
{
    unsigned error = lodepng::encode( pngFilename, image_,
                                      int( size_[0] ), int( size_[1] ),
                                      LCT_GREY, 8 );

    if ( error )
      throw std::invalid_argument( std::string( "Error while loading PNG file: " ) + lodepng_error_text(error)  );
}


double GrayScaleImage::operator() ( int x, int y ) const
{
    static const double maxVal = double( std::numeric_limits<unsigned char>::max() );
    return double ( getElement(x, y) ) / maxVal;
}

void GrayScaleImage::setElement( int x, int y, double val)
{
    assert( val <= 1.0 );
    assert( val >= 0.0 );
    getElement(x,y) = static_cast<unsigned char>( double( std::numeric_limits<unsigned char>::max() ) * val  );
}

GrayScaleImage GrayScaleImage::getResizedImage( unsigned int newWidth, unsigned int newHeight, bool bilinear ) const
{
    if ( newWidth == size_[0]  && newHeight == size_[1] )
      return *this;


    GrayScaleImage resizedImage;

    resizedImage.size_[0] = newWidth;
    resizedImage.size_[1] = newHeight;

    resizedImage.image_.resize( newWidth * newHeight );

    if ( bilinear )
    {
      double scaleX = double( size_[0]-1 ) / double( newWidth );
      double scaleY = double( size_[1]-1 ) / double( newHeight);

      for( int y = 0; y < int( newHeight ); ++y )
          for( int x = 0; x < int( newWidth ); ++x )
          {
            double oldX = double(x) * scaleX;
            double oldY = double(y) * scaleY;
            int oldXi = int( oldX );
            int oldYi = int( oldY );
            double xDiff = oldX - double(oldXi);
            double yDiff = oldY - double(oldYi);

            // bilinear interpolation

            resizedImage.getElement( x, y ) =
                      unsigned(
                      (1 - xDiff) * (1 - yDiff ) * getElement( oldXi    , oldYi    ) +
                          xDiff  * (1 - yDiff ) * getElement( oldXi + 1, oldYi    ) +
                      (1 - xDiff) *      yDiff   * getElement( oldXi    , oldYi + 1) +
                          xDiff  *      yDiff   * getElement( oldXi + 1, oldYi + 1) );
          }
    }
    else
    {
      double scaleX = double( size_[0] ) / double( newWidth );
      double scaleY = double( size_[1] ) / double( newHeight);

      for( int y = 0; y < int( newHeight ); ++y )
          for( int x = 0; x < int( newWidth ); ++x )
          {
            double oldX = double(x) * scaleX;
            double oldY = double(y) * scaleY;
            int oldXi = int( oldX );
            int oldYi = int( oldY );

            resizedImage.getElement( x, y ) = getElement( oldXi, oldYi );
          }
    }

    return resizedImage;
}

GrayScaleImage::pixel_t GrayScaleImage::pixelValueFromString( const std::string & str )
{
    std::stringstream ss (str);
    int tmp;
    ss >> std::hex >> tmp;
    return GrayScaleImage::pixel_t( tmp );
}








