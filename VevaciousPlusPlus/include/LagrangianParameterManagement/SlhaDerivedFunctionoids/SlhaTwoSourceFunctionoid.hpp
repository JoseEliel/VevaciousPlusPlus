/*
 * SlhaTwoSourceFunctionoid.hpp
 *
 *  Created on: Oct 27, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHATWOSOURCEFUNCTIONOID_HPP_
#define SLHATWOSOURCEFUNCTIONOID_HPP_

#include "CommonIncludes.hpp"
#include "LagrangianParameterManagement/SlhaDerivedFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class SlhaTwoSourceFunctionoid : public SlhaDerivedFunctionoid
  {
  public:
    SlhaTwoSourceFunctionoid(
            SlhaInterpolatedParameterFunctionoid const& firstChoiceFunctionoid,
         SlhaInterpolatedParameterFunctionoid const& secondChoiceFunctionoid );
    virtual ~SlhaTwoSourceFunctionoid();


    // This evaluates firstChoiceFunctionoid( logarithmOfScale ) and if it is
    // non-zero returns that value. Otherwise it returns
    // secondChoiceFunctionoid( logarithmOfScale ).
    virtual double operator()( double const logarithmOfScale ) const;

    // This returns the parameter value at firstChoiceIndex if it is non-zero,
    // otherwise it returns the parameter value at secondChoiceIndex.
    virtual double operator()( double const logarithmOfScale,
                       std::vector< double > const& interpolatedValues ) const
    { return ( ( interpolatedValues[ firstChoiceIndex ] == 0.0 ) ?
                   interpolatedValues[ secondChoiceIndex ] :
                   interpolatedValues[ firstChoiceIndex ] ); }

    // This is for creating a Python version of the potential.
    virtual std::string
    PythonParameterEvaluation( int const indentationSpaces ) const;


  protected:
    SlhaInterpolatedParameterFunctionoid const& firstChoiceFunctionoid;
    size_t const firstChoiceIndex;
    SlhaInterpolatedParameterFunctionoid const& secondChoiceFunctionoid;
    size_t const secondChoiceIndex;
  };





  // This evaluates firstChoiceFunctionoid( logarithmOfScale ) and if it is
  // non-zero returns that value. Otherwise it returns
  // secondChoiceFunctionoid( logarithmOfScale ).
  inline double
  SlhaTwoSourceFunctionoid::operator()( double const logarithmOfScale ) const
  {
    double firstChoiceValue( firstChoiceFunctionoid( logarithmOfScale ) );
    if( firstChoiceValue == 0.0 )
    {
      return secondChoiceFunctionoid( logarithmOfScale );
    }
    else
    {
      return firstChoiceValue;
    }
  }

  // This is for creating a Python version of the potential.
  inline std::string SlhaTwoSourceFunctionoid::PythonParameterEvaluation(
                                            int const indentationSpaces ) const
  {
    std::stringstream stringBuilder;
    stringBuilder << PythonIndent( indentationSpaces )
    << "parameterValues[ " << IndexInValuesVector()
    << " ] = FirstIfNonzeroOtherwiseSecond( parameterValues[ "
    << firstChoiceIndex
    << " ], parameterValues[ " << secondChoiceIndex << " ] )";
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */

#endif /* SLHATWOSOURCEFUNCTIONOID_HPP_ */
