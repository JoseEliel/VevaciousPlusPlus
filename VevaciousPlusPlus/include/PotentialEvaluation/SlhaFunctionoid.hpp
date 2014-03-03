/*
 * SlhaFunctionoid.hpp
 *
 *  Created on: Mar 3, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHAFUNCTIONOID_HPP_
#define SLHAFUNCTIONOID_HPP_

#include "../StandardIncludes.hpp"
#include "ParameterFunctionoid.hpp"
#include "SLHA.hpp"

namespace VevaciousPlusPlus
{

  class SlhaFunctionoid : public ParameterFunctionoid
  {
  public:
    SlhaFunctionoid( std::string const& indexString );
    virtual
    ~SlhaFunctionoid();


    void SetBlockPointer(
                LHPC::SLHA::SparseManyIndexedBlock< double >* slhaBlock )
    { this->slhaBlock = slhaBlock; }

    unsigned int NumberOfIndices() const{ return indexVector.size(); }

    // This updates currentValue based on logarithmOfScale.
    virtual void
    UpdateForNewLogarithmOfScale( double const logarithmOfScale );

    // This re-calculates the coefficients of the polynomial of the logarithm
    // of the scale used in evaluating the functionoid.
    void UpdateForNewSlhaParameters();

    // This is mainly for debugging.
    virtual std::string AsString();


  protected:
    LHPC::SLHA::SparseManyIndexedBlock< double >* slhaBlock;
    std::vector< int > const indexVector;
    std::vector< double > scaleLogarithmPowerCoefficients;
  };




  // This updates currentValue based on logarithmOfScale.
  inline void SlhaFunctionoid::UpdateForNewLogarithmOfScale(
                                                double const logarithmOfScale )
  {
    currentValue = scaleLogarithmPowerCoefficients[ 0 ];
    for( unsigned int whichPower( 1 );
         whichPower < scaleLogarithmPowerCoefficients.size();
         ++whichPower )
    {
      currentValue *= scaleLogarithmPowerCoefficients[ whichPower ];
      for( unsigned int powerCount( 0 );
           powerCount < whichPower;
           ++powerCount )
      {
        currentValue *= logarithmOfScale;
      }
    }

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "Block " << slhaBlock->getName() << " (for [";
    for( std::vector< int >::const_iterator
         whichIndex( indexVector.begin() );
         whichIndex < indexVector.end();
         ++whichIndex )
    {
      std::cout << " " << *whichIndex;
    }
    std::cout << " ])" << std::endl
    << "SlhaFunctionoid::UpdateForNewLogarithmOfScale( " << logarithmOfScale
    << " ) set currentValue to be " << currentValue;
    std::cout << std::endl;/**/
  }

  // This re-calculates the coefficients of the polynomial of the logarithm
  // of the scale used in evaluating the functionoid.
  inline void SlhaFunctionoid::UpdateForNewSlhaParameters()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "SlhaFunctionoid::UpdateForNewSlhaParameters()";
    std::cout << std::endl;/**/

    scaleLogarithmPowerCoefficients[ 0 ] = (*slhaBlock)[ 0 ]( indexVector );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "Block " << slhaBlock->getName() << " (for [";
    for( std::vector< int >::const_iterator
         whichIndex( indexVector.begin() );
         whichIndex < indexVector.end();
         ++whichIndex )
    {
      std::cout << ' ' << *whichIndex;
    }
    std::cout << " ])" << std::endl;
    for( int whichScale( 1 );
         whichScale <= slhaBlock->getNumberOfCopiesWithDifferentScale();
         ++whichScale )
    {
      std::cout
      << "scale copy " << whichScale << " has scale "
      << (*slhaBlock)[ whichScale ].getScale() << std::endl
      << "as string =" << std::endl
      << (*slhaBlock)[ whichScale ].getAsString();
    }
    std::cout << std::endl;/**/
  }

  // This is mainly for debugging.
  inline std::string SlhaFunctionoid::AsString()
  {
    std::stringstream returnStream;
    returnStream << slhaBlock->getName() << '[';
    for( std::vector< int >::const_iterator
         whichIndex( indexVector.begin() );
         whichIndex < indexVector.end();
         ++whichIndex )
    {
      returnStream << ' ' << *whichIndex;
    }
    returnStream << " ]";
    return std::string( returnStream.str() );
  }

} /* namespace VevaciousPlusPlus */
#endif /* SLHAFUNCTIONOID_HPP_ */
