/*
 * VectorUtilities.hpp
 *
 *  Created on: Nov 20, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef VECTORUTILITIES_HPP_
#define VECTORUTILITIES_HPP_

#include <vector>
#include <cstddef>
#include <cmath>

namespace VevaciousPlusPlus
{
  // This class is for a few little things with vectors that just crop up again
  // and again.
  class VectorUtilities
  {
  public:
    // This returns the Euclidean length squared.
    static double LengthSquared( std::vector< double > const& givenVector );

    // This returns true if each element in firstVector is within hypercubeSide
    // of the corresponding element in secondVector.
    static bool
    DifferenceIsWithinHypercube( std::vector< double > const& firstVector,
                                 std::vector< double > const& secondVector,
                                 double const hypercubeSide );
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

  // This returns true if each element in firstVector is within hypercubeSide
  // of the corresponding element in secondVector.
  inline bool VectorUtilities::DifferenceIsWithinHypercube(
                                      std::vector< double > const& firstVector,
                                     std::vector< double > const& secondVector,
                                                   double const hypercubeSide )
  {
    size_t const numberOfElements( firstVector.size() );
    if( secondVector.size() != numberOfElements )
    {
      return false;
    }
    for( size_t elementIndex( 0 );
         elementIndex < numberOfElements;
         ++elementIndex )
    {
      if( fabs( firstVector[ elementIndex ] - secondVector[ elementIndex ] )
          > hypercubeSide )
      {
        return false;
      }
    }
    return true;
  }

} /* namespace VevaciousPlusPlus */

#endif /* LHPC_VECTORUTILITIES_HPP_ */
