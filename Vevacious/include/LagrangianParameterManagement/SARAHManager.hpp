/*
 * SARAHManager.hpp
 *
 *  Created on: 9.4.2018
 *      Author: Simon Geisler 
 */

#ifndef SARAHMANAGER_HPP_
#define SARAHMANAGER_HPP_

#include "LesHouchesAccordBlockEntryManager.hpp"
#include <string>
#include <set>
#include <utility>
#include <cstddef>
#include <map>
#include <cstddef>
#include <utility>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <cmath>
#include <vector>
#include "LhaDerivedFunctionoids/LhaTwoSourceFunctionoid.hpp"
#include "LhaSourcedParameterFunctionoid.hpp"


namespace VevaciousPlusPlus
{

  class SARAHManager : public LesHouchesAccordBlockEntryManager
  {
  public:
    SARAHManager( std::string const& validBlocksString,
                                    std::string const& minimumScaleType,
                                    std::string const& minimumScaleArgument,
                                    std::string const& fixedScaleType,
                                    std::string const& fixedScaleArgument,
                                    std::string const& maximumScaleType,
                                    std::string const& maximumScaleArgument );
    SARAHManager(
                                 std::set< std::string > const& validBlocksSet,
                                    std::string const& minimumScaleType,
                                    std::string const& minimumScaleArgument,
                                    std::string const& fixedScaleType,
                                    std::string const& fixedScaleArgument,
                                    std::string const& maximumScaleType,
                                    std::string const& maximumScaleArgument );
    SARAHManager( std::string const& xmlFileName );
    virtual ~SARAHManager();


        // This fills the given vector with the values of the Lagrangian parameters
    // in activeInterpolatedParameters evaluated at the given scale, ordered so
    // that the indices given out by RegisterParameter correctly match the
    // parameter with its element in the vector. First it calls the base
    // version from LesHouchesAccordBlockEntryManager, then passes in the
    // vector in order to the functionoids in activeDerivedParameters, as each
    // derived parameter functionoid relies on the parameter values with
    // indices less than its own.
    virtual void ParameterValues( double const logarithmOfScale,
                              std::vector< double >& destinationVector ) const;
    // This first writes a function used by some derived parameters, and then
    // writes a function in the form
    // def LagrangianParameters( lnQ ): return ...
    // to return an array of the values of the Lagrangian parameters evaluated
    // at the scale exp(lnQ) (i.e. the logarithm of the scale is given as the
    // argument), in the order in which a call to ParametersAtScale would
    // return them internal to this C++ code.
    virtual std::string ParametersAsPython() const;

  protected:
      
    std::vector< LhaSourcedParameterFunctionoid* > activeDerivedParameters;
    std::map< std::string, std::string > aliasesToCaseStrings;
    
      // This adds newParameter to activeDerivedParameters and updates
    // numberOfDistinctActiveParameters and activeParametersToIndices (all
    // aliases in aliasesToCaseStrings which map to caseString are added
    // mapping to the index of this parameter), and returns true paired with
    // the index of newParameter.
    size_t AddNewDerivedParameter( std::string const& caseString,
                                LhaSourcedParameterFunctionoid* newParameter );

    // This pairs an index with true for convenience.
    std::pair< bool, size_t >
    PairAddNewDerivedParameter( std::string const& caseString,
                                LhaSourcedParameterFunctionoid* newParameter )
    { return std::pair< bool, size_t >( true,
                                        AddNewDerivedParameter( caseString,
                                                            newParameter ) ); }
     // This returns a string which is the concatenated set of strings from
    // parameter functionoids giving their Python evaluations.
    virtual std::string
    ParametersInPythonFunction( unsigned int const indentationSpaces ) const;
    //Add new Derived Parameters from Vector. Allows the use of IFNONZERO
    virtual void RegisterDerivedParameters(std::vector<std::pair<std::string,std::string>> derivedparameters);
  };
  // This fills the given vector with the values of the Lagrangian parameters
  // in activeInterpolatedParameters evaluated at the given scale, ordered so
  // that the indices given out by RegisterParameter correctly match the
  // parameter with its element in the vector. First it calls the base
  // version from LesHouchesAccordBlockEntryManager, then passes in the
  // vector in order to the functionoids in activeDerivedParameters, as each
  // derived parameter functionoid relies on the parameter values with
  // indices less than its own.
  inline void SARAHManager::ParameterValues(
                                                 double const logarithmOfScale,
                              std::vector< double >& destinationVector ) const
  {
    LesHouchesAccordBlockEntryManager::ParameterValues( logarithmOfScale,
                                                        destinationVector );
    for( std::vector< LhaSourcedParameterFunctionoid* >::const_iterator
         parameterInterpolator( activeDerivedParameters.begin() );
         parameterInterpolator < activeDerivedParameters.end();
         ++parameterInterpolator )
    {
      destinationVector[ (*parameterInterpolator)->IndexInValuesVector() ]
      = (*(*parameterInterpolator))( logarithmOfScale,
                                     destinationVector );
    }
  }

  // This returns a string which is the concatenated set of strings from
  // parameter functionoids giving their Python evaluations.
  inline std::string
  SARAHManager::ParametersInPythonFunction(
                                   unsigned int const indentationSpaces ) const
  {
    std::stringstream stringBuilder;
    stringBuilder
    << LesHouchesAccordBlockEntryManager::ParametersInPythonFunction(
                                                           indentationSpaces );
    for( std::vector< LhaSourcedParameterFunctionoid* >::const_iterator
         activeParameter( activeDerivedParameters.begin() );
         activeParameter < activeDerivedParameters.end();
         ++activeParameter )
    {
      stringBuilder
      << (*activeParameter)->PythonParameterEvaluation( indentationSpaces )
      << "\n";
    }
    return stringBuilder.str();
  }
  
  // This first writes a function used by some derived parameters, and then
  // writes a function in the form
  // def LagrangianParameters( lnQ ): return ...
  // to return an array of the values of the Lagrangian parameters evaluated
  // at the scale exp(lnQ) (i.e. the logarithm of the scale is given as the
  // argument), in the order in which a call to ParametersAtScale would
  // return them internal to this C++ code.
  inline std::string SARAHManager::ParametersAsPython() const
  {
    std::stringstream stringBuilder;
    stringBuilder
    << "def FirstIfNonzeroOtherwiseSecond( firstValue, secondValue ):\n"
    << "  if ( firstValue != 0.0 ):\n"
    << "    return firstValue\n"
    << "  else:\n"
    << "    return secondValue\n"
    << "\n"
    << LesHouchesAccordBlockEntryManager::ParametersAsPython();
    return stringBuilder.str();
  }
  
   inline size_t
  SARAHManager::AddNewDerivedParameter(
                                                 std::string const& caseString,
                                LhaSourcedParameterFunctionoid* newParameter )
  {
    activeDerivedParameters.push_back( newParameter );
    for( std::map< std::string, std::string >::const_iterator
         aliasToSwitchString( aliasesToCaseStrings.begin() );
         aliasToSwitchString != aliasesToCaseStrings.end();
         ++aliasToSwitchString )
    {
      if( aliasToSwitchString->second == caseString )
      {
        activeParametersToIndices[ aliasToSwitchString->first ]
        = newParameter->IndexInValuesVector();
      }
    }
    return RegisterNewParameter( *newParameter,
                                 caseString );
  }
} /* namespace VevaciousPlusPlus */

#endif /* SARAHMANAGER_HPP_ */
