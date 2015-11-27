/*
 * LhaTwoSourceFunctionoid.hpp
 *
 *  Created on: Oct 27, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef LHATWOSOURCEFUNCTIONOID_HPP_
#define LHATWOSOURCEFUNCTIONOID_HPP_

#include "CommonIncludes.hpp"
#include "LagrangianParameterManagement/LhaSourcedParameterFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class LhaTwoSourceFunctionoid : public LhaSourcedParameterFunctionoid
  {
  public:
    LhaTwoSourceFunctionoid( size_t const indexInValuesVector,
                              size_t const firstChoiceIndex,
                              size_t const secondChoiceIndex ) :
      LhaSourcedParameterFunctionoid( indexInValuesVector ),
      firstChoiceIndex( firstChoiceIndex ),
      secondChoiceIndex( secondChoiceIndex ) {}
    virtual ~LhaTwoSourceFunctionoid() {}


    // This returns the parameter value at firstChoiceIndex if it is non-zero,
    // otherwise it returns the parameter value at secondChoiceIndex.
    virtual double operator()( double const logarithmOfScale,
                        std::vector< double > const& interpolatedValues ) const
    { return ( ( interpolatedValues[ firstChoiceIndex ] != 0.0 ) ?
               interpolatedValues[ firstChoiceIndex ] :
               interpolatedValues[ secondChoiceIndex ] ); }

    // This returns the first parameter value if it is non-zero, otherwise the
    // second parameter value.
    double operator()( double const firstChoiceValue,
                       double const secondChoiceValue ) const
    { return
      ( ( firstChoiceValue != 0.0 ) ? firstChoiceValue : secondChoiceValue ); }

    // This is for creating a Python version of the potential.
    virtual std::string
    PythonParameterEvaluation( int const indentationSpaces ) const;

    // This is mainly for debugging.
    virtual std::string AsDebuggingString() const;


  protected:
    size_t const firstChoiceIndex;
    size_t const secondChoiceIndex;
  };





  // This is for creating a Python version of the potential.
  inline std::string LhaTwoSourceFunctionoid::PythonParameterEvaluation(
                                            int const indentationSpaces ) const
  {
    std::stringstream stringBuilder;
    stringBuilder << std::setprecision( 12 )
    << PythonIndent( indentationSpaces )
    << "parameterValues[ " << IndexInValuesVector()
    << " ] = FirstIfNonzeroOtherwiseSecond( parameterValues[ "
    << firstChoiceIndex
    << " ], parameterValues[ " << secondChoiceIndex << " ] )";
    return stringBuilder.str();
  }

  // This is mainly for debugging.
  inline std::string LhaTwoSourceFunctionoid::AsDebuggingString() const
  {
    std::stringstream stringBuilder;
    stringBuilder
    << "IndexInValuesVector() = " << IndexInValuesVector() << std::endl;
    stringBuilder << "firstChoiceIndex = " << firstChoiceIndex << std::endl;
    stringBuilder << "secondChoiceIndex = " << secondChoiceIndex << std::endl;
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */

#endif /* LHATWOSOURCEFUNCTIONOID_HPP_ */
