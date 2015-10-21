/*
 * LagrangianParameterManager.hpp
 *
 *  Created on: Oct 9, 2015
 *      Author: bol
 */

#ifndef LAGRANGIANPARAMETERMANAGER_HPP_
#define LAGRANGIANPARAMETERMANAGER_HPP_

#include "CommonIncludes.hpp"

namespace VevaciousPlusPlus
{

  class LagrangianParameterManager : public BOL::BasicObserved
  {
  public:
    LagrangianParameterManager();
    virtual ~LagrangianParameterManager();


    // This should ensure that the given string is understood as a Lagrangian
    // parameter and return an index which will be used to find the value of
    // the parameter in the vector returned by ParametersAtScale and in the
    // array returned by the function created by ParametersAsPython.
    virtual size_t RegisterParameter( std::string const& parameterName ) = 0;

    // This should return a vector of the values of the Lagrangian parameters
    // evaluated at the given scale, ordered so that the indices given out by
    // RegisterParameter correctly match the parameter with its element in the
    // returned vector.
    virtual std::vector< double >
    ParametersAtScale( double evaluationScale ) const = 0;

    // This puts all variables with index brackets into a consistent form. It
    // can be overridden if necessary.
    virtual std::string
    FormatVariable( std::string const& unformattedVariable ) const;

    // This should write a function in the form
    // def LagrangianParameter( Q ): return ...
    // to return an array of the values of the Lagrangian parameters evaluated
    // at the scale Q, in the order in which a call to ParametersAtScale would
    // return them internal to this C++ code.
    virtual std::string ParametersAsPython() const = 0;


  protected:
    // This puts all index brackets into a consistent form. It can be
    // overridden if necessary.
    virtual std::string FormatIndexBracketContent(
                          std::string const& unformattedBracketContent ) const;
  };





  // This puts all variables with index brackets into a consistent form.
  inline std::string LagrangianParameterManager::FormatVariable(
                                 std::string const& unformattedVariable ) const
  {
    size_t openBracket( unformattedVariable.find( '[' ) );
    if( openBracket == std::string::npos )
    {
      return std::string( unformattedVariable );
    }
    if( unformattedVariable[ unformattedVariable.size() - 1 ] != ']' )
    {
      throw std::runtime_error(
                         "In parsing model file, [...] not closed properly." );
    }
    std::stringstream formattedStream;
    formattedStream << unformattedVariable.substr( 0,
                                                   openBracket ) << '['
    << FormatIndexBracketContent( unformattedVariable.substr(
                                                           ( openBracket + 1 ),
                   ( unformattedVariable.size() - openBracket - 2 ) ) ) << ']';
    return formattedStream.str();
  }

  // This puts all index brackets into a consistent form.
  inline std::string LagrangianParameterManager::FormatIndexBracketContent(
                           std::string const& unformattedBracketContent ) const
  {
    std::vector< int > indicesVector( BOL::StringParser::stringToIntVector(
                                                 unformattedBracketContent ) );
    std::stringstream indicesStream;
    for( std::vector< int >::iterator
         whichIndex( indicesVector.begin() );
         whichIndex < indicesVector.end();
         ++whichIndex )
    {
      if( whichIndex != indicesVector.begin() )
      {
        indicesStream << ',';
      }
      indicesStream << *whichIndex;
    }
    return indicesStream.str();
  }

} /* namespace VevaciousPlusPlus */

#endif /* LAGRANGIANPARAMETERMANAGER_HPP_ */
