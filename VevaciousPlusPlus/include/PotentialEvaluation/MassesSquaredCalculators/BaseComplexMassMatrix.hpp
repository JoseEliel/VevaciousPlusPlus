/*
 * BaseComplexMassMatrix.hpp
 *
 *  Created on: Oct 09, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BASECOMPLEXMASSMATRIX_HPP_
#define BASECOMPLEXMASSMATRIX_HPP_

#include "MassesSquaredFromMatrix.hpp"
#include <complex>
#include "PotentialEvaluation/BuildingBlocks/ParametersAndFieldsProductSum.hpp"
#include <map>
#include <string>
#include <vector>
#include <sstream>

namespace VevaciousPlusPlus
{

  class BaseComplexMassMatrix :
                       public MassesSquaredFromMatrix< std::complex< double > >
  {
  public:
    typedef
    std::pair< ParametersAndFieldsProductSum, ParametersAndFieldsProductSum >
    ComplexParametersAndFieldsProductSum;

    BaseComplexMassMatrix( size_t const numberOfElements,
                   std::map< std::string, std::string > const& attributeMap ) :
      MassesSquaredFromMatrix< std::complex< double > >( numberOfRows,
                                                         attributeMap ),
      matrixElements( ( numberOfRows * numberOfRows ),
         ComplexParametersAndFieldsProductSum( ParametersAndFieldsProductSum(),
                                         ParametersAndFieldsProductSum() ) ) {}

    BaseComplexMassMatrix( BaseComplexMassMatrix const& copySource ) :
      MassesSquaredFromMatrix< std::complex< double > >( copySource ),
      matrixElements( copySource.matrixElements ) {}

    BaseComplexMassMatrix() :
      MassesSquaredFromMatrix< std::complex< double > >(),
      matrixElements() {}

    virtual ~BaseComplexMassMatrix() {}


    // This calls UpdateForFixedScale on each element of matrixElements.
    virtual void
    UpdateForFixedScale( std::vector< double > const& parameterValues );

    // This allows access to the pair of polynomial sums for a given index.
    ComplexParametersAndFieldsProductSum&
    ElementAt( size_t const elementIndex )
    { return matrixElements[ elementIndex ]; }

    // This is mainly for debugging:
    std::string AsString() const;

    std::vector< ComplexParametersAndFieldsProductSum > const&
    MatrixElements() const{ return matrixElements; }


  protected:
    std::vector< ComplexParametersAndFieldsProductSum > matrixElements;
  };





  // This calls UpdateForFixedScale on each element of matrixElements.
  inline void BaseComplexMassMatrix::UpdateForFixedScale(
                                 std::vector< double > const& parameterValues )
  {
    for( std::vector< ComplexParametersAndFieldsProductSum >::iterator
         complexPair( matrixElements.begin() );
         complexPair < matrixElements.end();
         ++complexPair )
    {
      complexPair->first.UpdateForFixedScale( parameterValues );
      complexPair->second.UpdateForFixedScale( parameterValues );
    }
  }

  // This is mainly for debugging:
  inline std::string BaseComplexMassMatrix::AsString() const
  {
    std::stringstream returnStream;
    returnStream
    << "numberOfRows = " << numberOfRows << ", matrixElements = ";
    for( std::vector< ComplexParametersAndFieldsProductSum >::const_iterator
         whichPair( matrixElements.begin() );
         whichPair < matrixElements.end();
         ++whichPair )
    {
      returnStream
      << std::endl << "real:" << std::endl
      << whichPair->first.AsDebuggingString()
      << std::endl << "imag:" << std::endl
      << whichPair->second.AsDebuggingString()
      << std::endl;
    }
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */

#endif /* BASECOMPLEXMASSMATRIX_HPP_ */
