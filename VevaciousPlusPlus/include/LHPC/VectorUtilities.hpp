/*
 * VectorUtilities.hpp
 *
 *  Created on: Nov 20, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef LHPC_VECTORUTILITIES_HPP_
#define LHPC_VECTORUTILITIES_HPP_

#include <vector>
#include <sstream>

namespace LHPC
{
  // This class is for a few little things with vectors that just crop up again
  // and again.
  class VectorUtilities
  {
  public:
    // This returns the Euclidean length squared.
    static double LengthSquared( std::vector< double > const& givenVector );
  };





  // This returns the Euclidean length squared.
  inline double VectorUtilities::LengthSquared(
                                     std::vector< double > const& givenVector )
  {
    double lengthSquared( 0.0 );
    for( std::vector< double >::const_iterator
         vectorElement( givenVector.begin() );
         vectorElement != givenVector.end();
         ++vectorElement )
    {
      lengthSquared += ( (*vectorElement) * (*vectorElement) );
    }
    return lengthSquared;
  }
}

#endif /* LHPC_VECTORUTILITIES_HPP_ */
