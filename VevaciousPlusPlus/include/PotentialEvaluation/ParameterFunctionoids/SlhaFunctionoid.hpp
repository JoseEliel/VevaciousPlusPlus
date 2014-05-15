/*
 * SlhaFunctionoid.hpp
 *
 *  Created on: Mar 3, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHAFUNCTIONOID_HPP_
#define SLHAFUNCTIONOID_HPP_

#include "../../CommonIncludes.hpp"
#include "ParameterFunctionoid.hpp"
#include "Eigen/Dense"
#include "../SimplePolynomial.hpp"


namespace VevaciousPlusPlus
{

  class SlhaFunctionoid : public ParameterFunctionoid
  {
  public:
    SlhaFunctionoid( std::string const& indexString,
                     std::string const& creationString,
                     std::string const& pythonParameterName );
    virtual
    ~SlhaFunctionoid();


    void SetBlockPointer(
                LHPC::SLHA::SparseManyIndexedBlock< double >* slhaBlock )
    { this->slhaBlock = slhaBlock; }

    unsigned int NumberOfIndices() const{ return indexVector.size(); }

    // This returns the value of the functionoid for the given logarithm of the
    // scale.
    virtual double operator()( double const logarithmOfScale ) const
    { return scaleLogarithmPowerCoefficients( logarithmOfScale ); }

    // This updates currentValue based on logarithmOfScale.
    virtual void
    UpdateForNewLogarithmOfScale( double const logarithmOfScale )
    { currentValue = (*this)( logarithmOfScale ); }

    // This re-calculates the coefficients of the polynomial of the logarithm
    // of the scale used in evaluating the functionoid.
    void UpdateForNewSlhaParameters();

    // This is mainly for debugging.
    virtual std::string AsString();

    // This is for creating a Python version of the potential.
    virtual std::string PythonParameterEvaluation() const;


  protected:
    LHPC::SLHA::SparseManyIndexedBlock< double >* slhaBlock;
    std::vector< int > const indexVector;
    SimplePolynomial scaleLogarithmPowerCoefficients;
  };




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
    return returnStream.str();
  }

  // This is for creating a Python version of the potential.
  inline std::string SlhaFunctionoid::PythonParameterEvaluation() const
  {
    std::vector< double > const&
    scaleCoefficients( scaleLogarithmPowerCoefficients.CoefficientVector() );
    std::stringstream stringBuilder;
    stringBuilder << std::setprecision( 12 ) << pythonParameterName << " = ( "
    << scaleCoefficients[ 0 ];
    for( unsigned int whichPower( 1 );
         whichPower < scaleCoefficients.size();
         ++whichPower )
    {
      stringBuilder << " + ( " << scaleCoefficients[ whichPower ]
      << " ) * lnQ**" << whichPower;
    }
    stringBuilder << " )";
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* SLHAFUNCTIONOID_HPP_ */
