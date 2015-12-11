/*
 * MinuitHypersphereBoundAlternative.hpp
 *
 *  Created on: Oct 6, 2015
 *      Author: bol
 */

#ifndef MINUITHYPERSPHEREBOUNDALTERNATIVE_HPP_
#define MINUITHYPERSPHEREBOUNDALTERNATIVE_HPP_

#include <vector>
#include <cmath>

namespace VevaciousPlusPlus
{

  class MinuitHypersphereBoundAlternative
  {
  public:
    // This scales down the values of variableVector if the sum of the squares
    // of the elements is larger than capLengthSquared so that the sum is equal
    // to capLengthSquared, and returns the difference of the original sum of
    // squares from capLengthSquared.
    static double CapVariableVector( std::vector< double >& variableVector,
                                     double const capLengthSquared );
  };




  // This scales down the values of variableVector if the sum of the squares of
  // the elements is larger than capLengthSquared so that the sum is equal to
  // capLengthSquared, and returns the difference of the original sum of
  // squares from capLengthSquared.
  inline double MinuitHypersphereBoundAlternative::CapVariableVector(
                                         std::vector< double >& variableVector,
                                                double const capLengthSquared )
  {
    double lengthSquared( 0.0 );
    for( std::vector< double >::const_iterator
         variableValue( variableVector.begin() );
         variableValue < variableVector.end();
         ++variableValue )
    {
      lengthSquared += ( (*variableValue) * (*variableValue) );
    }
    if( lengthSquared <= capLengthSquared )
    {
      return 0.0;
    }
    double const scaleFactor( sqrt( capLengthSquared / lengthSquared ) );
    for( size_t fieldIndex( 0 );
         fieldIndex < variableVector.size();
         ++fieldIndex )
    {
      variableVector[ fieldIndex ] *= scaleFactor;
    }
    return ( lengthSquared - capLengthSquared );
  }

} /* namespace VevaciousPlusPlus */

#endif /* MINUITHYPERSPHEREBOUNDALTERNATIVE_HPP_ */
