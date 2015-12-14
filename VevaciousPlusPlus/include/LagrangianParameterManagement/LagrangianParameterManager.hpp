/*
 * LagrangianParameterManager.hpp
 *
 *  Created on: Oct 9, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef LAGRANGIANPARAMETERMANAGER_HPP_
#define LAGRANGIANPARAMETERMANAGER_HPP_

#include "LHPC/Utilities/BasicObserverPattern.hpp"
#include <string>
#include <utility>
#include <cstddef>
#include <vector>
#include <sstream>
#include "LHPC/Utilities/ParsingUtilities.hpp"
#include <stdexcept>

namespace VevaciousPlusPlus
{

  class LagrangianParameterManager : public LHPC::BasicObserved
  {
  public:
    LagrangianParameterManager() : LHPC::BasicObserved() {}
    virtual ~LagrangianParameterManager() {}


    // This should return the value of the requested parameter at the requested
    // scale (exp( logarithmOfScale )) as an alternative to adding it to the
    // list of parameters evaluated by ParameterValues. It is almost certain to
    // be much slower if used to obtain parameters repeatedly for different
    // scales at the same parameter point, but is more efficient for the rest
    // of the execution of Vevacious, which depends on the speed of evaluation
    // of the potential, if there are some parameters which do not need to be
    // evaluated for the potential but still depend on Lagrangian parameters,
    // for example when evaluating the VEVs of the DSB vacuum for the parameter
    // point.
    virtual double OnceOffParameter( std::string const& parameterName,
                                     double const logarithmOfScale ) const = 0;

    // This should check whether the given string is understood as a Lagrangian
    // parameter and return a pair of the bool indicating whether it was a
    // valid string for an understood Lagrangian parameter, paired with an
    // index which will be used to find the value of the parameter in the
    // vector returned by ParametersAtScale and in the array returned by the
    // function created by ParametersAsPython. It would be inconvenient to
    // throw an exception for an invalid parameter input string as it is
    // expected that sometimes strings for fields will be passed in when
    // parsing a model file where certain fields have been de-selected as
    // variable fields. Hence calls of this method would always be wrapped in
    // try-catch, making the whole thing a bit more complicated than it has to
    // be, in my opinion.
    virtual std::pair< bool, size_t >
    RegisterParameter( std::string const& parameterName ) = 0;

    // This should fill the given vector with the values of the Lagrangian
    // parameters evaluated at the given scale, ordered so that the indices
    // given out by RegisterParameter correctly match the parameter with its
    // element in the vector.
    virtual void ParameterValues( double logarithmOfScale,
                          std::vector< double >& destinationVector ) const = 0;

    // This should return the minimum scale which is appropriate for evaluating
    // the Lagrangian parameters at the current parameter point.
    virtual double MinimumEvaluationScale() const = 0;

    // This should return a scale which is appropriate for using for a
    // fixed-scale calculation for the current parameter point.
    virtual double AppropriateSingleFixedScale() const = 0;

    // This should return the minimum scale which is appropriate for evaluating
    // the Lagrangian parameters at the current parameter point.
    virtual double MaximumEvaluationScale() const = 0;

    // This just runs the internal PrepareNewParameterPoint method then updates
    // the observers.
    virtual void NewParameterPoint( std::string const& newInput );

    // This puts all variables with index brackets into a consistent form. It
    // can be overridden if necessary.
    virtual std::string
    FormatVariable( std::string const& variableToFormat ) const;

    // This should write a function in the form
    // def LagrangianParameters( lnQ ): return ...
    // to return an array of the values of the Lagrangian parameters evaluated
    // at the scale exp(lnQ) (i.e. the logarithm of the scale is given as the
    // argument), in the order in which a call to ParametersAtScale would
    // return them internal to this C++ code. By default it produces Python
    // code which throws an exception, as a derived class may not need to write
    // Python as the CosmoTransitions Python code may not be needed.
    virtual std::string ParametersAsPython() const
    { return std::string( "def LagrangianParameters( lnQ ): raise"
                          " NotImplementedError( \"A C++ class derived from"
                          " Vevacious::LagrangainParameterManager did not"
                          " override the ParametersAsPython() method to"
                          " provide valid Python code to return the values of"
                          " the Lagrangian parameters at the scale given by"
                          " exp(lnQ) as an array.\")");}


  protected:
    // This should prepare the LagrangianParameterManager for a new parameter
    // point. It is expected that newInput is the name of a file containing
    // input required to evaluate the Lagrangian parameters (even the values
    // of the parameters themselves, such as given by an SLHA file), but
    // could itself be an encoding of the new input directly if a derived class
    // needs it. It could even just be ignored, if the derived class takes on
    // the role of deciding parameter input itself. Regardless, there should be
    // some code to prepare the parameters so that code objects which depend on
    // them can be informed that there is a new parameter point.
    virtual void PrepareNewParameterPoint( std::string const& newInput ) = 0;

    // This puts all index brackets into a consistent form. It can be
    // overridden if necessary.
    virtual std::string FormatIndexBracketContent(
                          std::string const& unformattedBracketContent ) const;

    // This returns "everything up to the '['" paired with "everything within
    // the '[' to the ']' (not including the brackets themselves), throwing
    // exceptions if there anything other than either no square brackets, or a
    // single '[' with a single ']' coming after it.
    std::pair< std::string, std::string >
    SeparateIndexBracket( std::string const& stringToSeparate ) const;
  };





  // This just runs the internal PrepareNewParameterPoint method then updates
  // the observers.
  inline void
  LagrangianParameterManager::NewParameterPoint( std::string const& newInput )
  {
    PrepareNewParameterPoint( newInput );
    UpdateObservers();
  }

  // This puts all variables with index brackets into a consistent form.
  inline std::string LagrangianParameterManager::FormatVariable(
                                    std::string const& variableToFormat ) const
  {
    std::pair< std::string, std::string > const
    separatedVariable( SeparateIndexBracket( variableToFormat ) );
    if( separatedVariable.second.empty() )
    {
      return separatedVariable.first;
    }
    std::stringstream formattedStream;
    formattedStream << separatedVariable.first << '['
    << FormatIndexBracketContent( separatedVariable.second ) << ']';
    return formattedStream.str();
  }

  // This puts all index brackets into a consistent form.
  inline std::string LagrangianParameterManager::FormatIndexBracketContent(
                           std::string const& unformattedBracketContent ) const
  {
    std::vector< std::string > const
    indicesVector( LHPC::ParsingUtilities::SplitBySubstrings(
                                                     unformattedBracketContent,
                                                              " \t\n\r,;" ) );
    std::stringstream indicesStream;
    for( std::vector< std::string >::const_iterator
         whichIndex( indicesVector.begin() );
         whichIndex < indicesVector.end();
         ++whichIndex )
    {
      if( whichIndex != indicesVector.begin() )
      {
        indicesStream << ',';
      }
      indicesStream << std::atoi( whichIndex->c_str() );
    }
    return indicesStream.str();
  }

  // This returns "everything up to the '['" paired with "everything within the
  // '[' to the ']' (not including the brackets themselves), throwing
  // exceptions if there anything other than either no square brackets, or a
  // single '[' with a single ']' coming after it.
  inline std::pair< std::string, std::string >
  LagrangianParameterManager::SeparateIndexBracket(
                                    std::string const& stringToSeparate ) const
  {
    size_t const openBracket( stringToSeparate.find( '[' ) );
    if( openBracket == std::string::npos )
    {
      return std::make_pair( stringToSeparate, std::string( "" ) );
    }
    size_t const closeBracket( stringToSeparate.find_last_of( '[' ) );
    if( ( closeBracket == std::string::npos )
        ||
        ( stringToSeparate.find( openBracket, '[' ) != std::string::npos )
        ||
        ( stringToSeparate.find_last_of( closeBracket,
                                         '[' ) != std::string::npos ) )
    {
      throw std::runtime_error(
                           "In parsing variable, [...] not closed properly." );
    }
    return std::make_pair( stringToSeparate.substr( 0, openBracket ),
                           stringToSeparate.substr( ( openBracket + 1 ),
                             ( stringToSeparate.size() - openBracket - 2 ) ) );
  }

} /* namespace VevaciousPlusPlus */

#endif /* LAGRANGIANPARAMETERMANAGER_HPP_ */
