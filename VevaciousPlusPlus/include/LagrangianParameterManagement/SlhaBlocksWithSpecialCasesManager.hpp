/*
 * SlhaBlocksWithSpecialCasesManager.hpp
 *
 *  Created on: Nov 2, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHABLOCKSWITHSPECIALCASESMANAGER_HPP_
#define SLHABLOCKSWITHSPECIALCASESMANAGER_HPP_

#include "CommonIncludes.hpp"
#include "LesHouchesAccordBlockEntryManager.hpp"
#include "LhaSourcedParameterFunctionoid.hpp"
#include "SlhaDerivedFunctionoids/SlhaDsbHiggsVevFunctionoid.hpp"
#include "SlhaDerivedFunctionoids/SlhaHiggsMixingBilinearFunctionoid.hpp"
#include "SlhaDerivedFunctionoids/SlhaMassSquaredDiagonalFunctionoid.hpp"
#include "SlhaDerivedFunctionoids/SlhaTrilinearDiagonalFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class SlhaBlocksWithSpecialCasesManager :
                                       public LesHouchesAccordBlockEntryManager
  {
  public:
    SlhaBlocksWithSpecialCasesManager( std::string const& validBlocksString,
                                       std::string const& minimumScaleType,
                                       std::string const& minimumScaleArgument,
                                       std::string const& fixedScaleType,
                                       std::string const& fixedScaleArgument,
                                       std::string const& maximumScaleType,
                                     std::string const& maximumScaleArgument );
    SlhaBlocksWithSpecialCasesManager(
                                 std::set< std::string > const& validBlocksSet,
                                       std::string const& minimumScaleType,
                                       std::string const& minimumScaleArgument,
                                       std::string const& fixedScaleType,
                                       std::string const& fixedScaleArgument,
                                       std::string const& maximumScaleType,
                                     std::string const& maximumScaleArgument );
    SlhaBlocksWithSpecialCasesManager( std::string const& xmlFileName );
    virtual ~SlhaBlocksWithSpecialCasesManager();


    // This returns the value of the requested parameter at the requested scale
    // (exp( logarithmOfScale )) as an alternative to adding it to the list of
    // parameters evaluated by ParameterValues. It is almost certain to
    // be much slower (as it searches for the appropriate block and then makes
    // a new functionoid) if used to obtain parameters repeatedly for
    // different scales at the same parameter point, but is more efficient for
    // the rest of the execution of Vevacious, which depends on the speed of
    // evaluation of the potential, if there are some parameters which do not
    // need to be evaluated for the potential but still depend on Lagrangian
    // parameters, for example when evaluating the VEVs of the DSB vacuum for
    // the parameter point.
    virtual double OnceOffParameter( std::string const& parameterName,
                                     double const logarithmOfScale ) const;

    // This fills the given vector with the values of the Lagrangian parameters
    // in activeInterpolatedParameters evaluated at the given scale, ordered so
    // that the indices given out by RegisterParameter correctly match the
    // parameter with its element in the vector. First it calls the base
    // version from LesHouchesAccordBlockEntryManager, then passes in the
    // vector in order to the functionoids in activeDerivedParameters, as each
    // derived parameter functionoid relies on the parameter values with
    // indices less than its own.
    virtual void
    ParameterValues( double const logarithmOfScale,
                     std::vector< double >& destinationVector ) const;

    // This first writes a function used by some derived parameters, and then
    // writes a function in the form
    // def LagrangianParameters( lnQ ): return ...
    // to return an array of the values of the Lagrangian parameters evaluated
    // at the scale exp(lnQ) (i.e. the logarithm of the scale is given as the
    // argument), in the order in which a call to ParametersAtScale would
    // return them internal to this C++ code.
    virtual std::string ParametersAsPython() const;

    // This is mainly for debugging.
    virtual std::string AsDebuggingString() const;


  protected:
    static bool
    SortParameterByIndex( LhaSourcedParameterFunctionoid const* firstPointer,
                         LhaSourcedParameterFunctionoid const* secondPointer )
    { return ( firstPointer->IndexInValuesVector()
               < secondPointer->IndexInValuesVector() ); }

    std::vector< LhaSourcedParameterFunctionoid* > activeDerivedParameters;
    std::map< std::string, std::string > aliasesToCaseStrings;

    // This adds all the valid aliases to aliasesToSwitchStrings.
    void InitializeSlhaOneOrTwoAliases();

    // This duplicates a lot of code from RegisterUnregisteredSpecialCase, but
    // there doesn't seem to be an elegant way of using the common code as
    // there is too much entanglement with registering new parameters or not.
    virtual double OnceOffSpecialCase( std::string const& parameterName,
                                       double const logarithmOfScale ) const;

    // This adds newParameter to activeDerivedParameters and updates
    // numberOfDistinctActiveParameters and activeParametersToIndices (all
    // aliases in aliasesToCaseStrings which map to caseString are added
    // mapping to the index of this parameter), and returns true paired with
    // the index of newParameter.
    std::pair< bool, size_t >
    AddNewDerivedParameter( std::string const& caseString,
                            LhaSourcedParameterFunctionoid* newParameter );

    // This checks parameterName against all the special cases. If the
    // parameter name is not recognized as a valid special case, then the
    // LesHouchesAccordBlockEntryManager result is returned, otherwise the
    // special case is registered and returned.
    virtual std::pair< bool, size_t >
    RegisterUnregisteredParameter( std::string const& parameterName );

    // This adds the parameter based on the alias given by switchString for the
    // parameter.
    virtual std::pair< bool, size_t >
    RegisterUnregisteredSpecialCase( std::string const& caseString );

    // This registers the required TX[0,0], AX[0,0], and YX[0,0] parameters
    // where X is replaced by sfermionType and 0 by generationIndex, and then
    // returns the result of AddNewDerivedParameter using a new
    // SlhaTrilinearDiagonalFunctionoid constructed with those parameter
    // functionoids.
    std::pair< bool, size_t >
    RegisterSlhaOneOrTwoCompatibleTrilinear( std::string const& caseString,
                                             char const sfermionType,
                                             char const generationIndex );

    // This creates a temporary functionoid and passes it parameter values from
    // OnceOff(...) to return the appropriate value for caseString.
    double
    OnceOffSlhaOneOrTwoCompatibleTrilinear( std::string const& caseString,
                                            char const sfermionType,
                                            char const generationIndex,
                                         double const logarithmOfScale ) const;

    // This registers the required M2X[0,0] and MSOFT[msoftIndex] parameters
    // where X is replaced by sfermionType and 0 by generationIndex, and then
    // returns the result of AddNewDerivedParameter using a new
    // SlhaMassSquaredDiagonalFunctionoid constructed with those parameter
    // functionoids.
    std::pair< bool, size_t >
    RegisterSlhaOneOrTwoCompatibleMassSquared( std::string const& caseString,
                                               char const sfermionType,
                                               char const generationIndex,
                                               unsigned int const msoftIndex );

    // This creates a temporary functionoid and passes it parameter values from
    // OnceOff(...) to return the appropriate value for caseString.
    double
    OnceOffSlhaOneOrTwoCompatibleMassSquared( std::string const& caseString,
                                              char const sfermionType,
                                              char const generationIndex,
                                              unsigned int const msoftIndex,
                                         double const logarithmOfScale ) const;

    // This maps both caseString and the format-corrected version of
    // blockString to caseString in aliasesToSwitchStrings.
    void MapCaseStringAndSlhaBlockToCaseString( std::string const& caseString,
                                               std::string const& blockString )
    { aliasesToCaseStrings[ caseString ]
      = aliasesToCaseStrings[ FormatVariable( blockString ) ] = caseString; }

    // This returns a string which is the concatenated set of strings from
    // parameter functionoids giving their Python evaluations.
    virtual std::string
    ParametersInPythonFunction( unsigned int const indentationSpaces ) const;
  };





  // This returns the value of the requested parameter at the requested scale
  // (exp( logarithmOfScale )) as an alternative to adding it to the list of
  // parameters evaluated by ParameterValues. It is almost certain to
  // be much slower (as it searches for the appropriate block and then makes
  // a new functionoid) if used to obtain parameters repeatedly for
  // different scales at the same parameter point, but is more efficient for
  // the rest of the execution of Vevacious, which depends on the speed of
  // evaluation of the potential, if there are some parameters which do not
  // need to be evaluated for the potential but still depend on Lagrangian
  // parameters, for example when evaluating the VEVs of the DSB vacuum for
  // the parameter point.
  inline double SlhaBlocksWithSpecialCasesManager::OnceOffParameter(
                                              std::string const& parameterName,
                                          double const logarithmOfScale ) const
  {
    std::map< std::string, std::string >::const_iterator
    aliasToSwitchString( aliasesToCaseStrings.find( parameterName ) );
    if( aliasToSwitchString != aliasesToCaseStrings.end() )
    {
      return OnceOffSpecialCase( aliasToSwitchString->second,
                                 logarithmOfScale );
    }
    else
    {
      return
      LesHouchesAccordBlockEntryManager::OnceOffParameter( parameterName,
                                                           logarithmOfScale );
    }
  }

  // This fills the given vector with the values of the Lagrangian parameters
  // in activeInterpolatedParameters evaluated at the given scale, ordered so
  // that the indices given out by RegisterParameter correctly match the
  // parameter with its element in the vector. First it calls the base
  // version from LesHouchesAccordBlockEntryManager, then passes in the
  // vector in order to the functionoids in activeDerivedParameters, as each
  // derived parameter functionoid relies on the parameter values with
  // indices less than its own.
  inline void SlhaBlocksWithSpecialCasesManager::ParameterValues(
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

  // This first writes a function used by some derived parameters, and then
  // writes a function in the form
  // def LagrangianParameters( lnQ ): return ...
  // to return an array of the values of the Lagrangian parameters evaluated
  // at the scale exp(lnQ) (i.e. the logarithm of the scale is given as the
  // argument), in the order in which a call to ParametersAtScale would
  // return them internal to this C++ code.
  inline std::string
  SlhaBlocksWithSpecialCasesManager::ParametersAsPython() const
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

  // This adds newParameter to activeDerivedParameters and updates
  // numberOfDistinctActiveParameters and activeParametersToIndices (all
  // aliases in aliasesToSwitchStrings which map to switchString are added
  // mapping to the index of this parameter), and returns true paired with
  // the index of newParameter.
  inline std::pair< bool, size_t >
  SlhaBlocksWithSpecialCasesManager::AddNewDerivedParameter(
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

  // This checks parameterName against all the special cases. If the
  // parameter name is not recognized as a valid special case, then the
  // LesHouchesAccordBlockEntryManager result is returned, otherwise the
  // special case is registered and returned.
  inline std::pair< bool, size_t >
  SlhaBlocksWithSpecialCasesManager::RegisterUnregisteredParameter(
                                             std::string const& parameterName )
  {
    std::map< std::string, std::string >::const_iterator
    aliasToSwitchString( aliasesToCaseStrings.find( parameterName ) );
    if( aliasToSwitchString != aliasesToCaseStrings.end() )
    {
      return RegisterUnregisteredSpecialCase( aliasToSwitchString->second );
    }
    else
    {
      return LesHouchesAccordBlockEntryManager::RegisterUnregisteredParameter(
                                                               parameterName );
    }
  }

  // This registers the required TX[0,0], AX[0,0], and YX[0,0] parameters
  // where X is replaced by sfermionType and 0 by generationIndex, and then
  // returns the result of AddNewDerivedParameter using a new
  // SlhaTrilinearDiagonalFunctionoid constructed with those parameter
  // functionoids.
  inline std::pair< bool, size_t >
  SlhaBlocksWithSpecialCasesManager::RegisterSlhaOneOrTwoCompatibleTrilinear(
                                                 std::string const& caseString,
                                                       char const sfermionType,
                                                   char const generationIndex )
  {
    std::string blockNameWithIndices( "T?[?,?]" );
    blockNameWithIndices[ 1 ] = sfermionType;
    blockNameWithIndices[ 3 ] = blockNameWithIndices[ 5 ] = generationIndex;
    size_t const
    directTrilinearIndex( RegisterBlockEntry( blockNameWithIndices ) );
    blockNameWithIndices[ 0 ] = 'A';
    size_t const
    trilinearOverYukawaIndex( RegisterBlockEntry( blockNameWithIndices ) );
    blockNameWithIndices[ 0 ] = 'Y';
    size_t const
    appropriateYukawaIndex( RegisterBlockEntry( blockNameWithIndices ) );
    return AddNewDerivedParameter( caseString,
                                   new SlhaTrilinearDiagonalFunctionoid(
                                              numberOfDistinctActiveParameters,
                                                          directTrilinearIndex,
                                                      trilinearOverYukawaIndex,
                                                    appropriateYukawaIndex ) );
  }

  // This creates a temporary functionoid and passes it parameter values from
  // OnceOff(...) to return the appropriate value for caseString.
  inline double
  SlhaBlocksWithSpecialCasesManager::OnceOffSlhaOneOrTwoCompatibleTrilinear(
                                                 std::string const& caseString,
                                                       char const sfermionType,
                                                    char const generationIndex,
                                          double const logarithmOfScale ) const
  {
    std::string blockNameWithIndices( "T?[?,?]" );
    blockNameWithIndices[ 1 ] = sfermionType;
    blockNameWithIndices[ 3 ] = blockNameWithIndices[ 5 ] = generationIndex;
    double const directTrilinear( OnceOffParameter( blockNameWithIndices,
                                                    logarithmOfScale ) );
    blockNameWithIndices[ 0 ] = 'A';
    double const trilinearOverYukawa( OnceOffParameter( blockNameWithIndices,
                                                    logarithmOfScale ) );
    blockNameWithIndices[ 0 ] = 'Y';
    double const appropriateYukawa( OnceOffParameter( blockNameWithIndices,
                                                      logarithmOfScale ) );

    SlhaTrilinearDiagonalFunctionoid temporaryParameter( 0,
                                                         0,
                                                         0,
                                                         0 );
    return temporaryParameter( directTrilinear,
                               trilinearOverYukawa,
                               appropriateYukawa );
  }

  // This registers the required M2X[0,0] and MSOFT[msoftIndex] parameters
  // where X is replaced by sfermionType and 0 by generationIndex, and then
  // returns the result of AddNewDerivedParameter using a new
  // SlhaMassSquaredDiagonalFunctionoid constructed with those parameter
  // functionoids.
  inline std::pair< bool, size_t >
  SlhaBlocksWithSpecialCasesManager::RegisterSlhaOneOrTwoCompatibleMassSquared(
                                                 std::string const& caseString,
                                                       char const sfermionType,
                                                    char const generationIndex,
                                                unsigned int const msoftIndex )
  {
    std::string blockNameWithIndices( "MS?2[?,?]" );
    blockNameWithIndices[ 2 ] = sfermionType;
    blockNameWithIndices[ 5 ] = blockNameWithIndices[ 7 ] = generationIndex;
    size_t const squareMassIndex( RegisterBlockEntry(blockNameWithIndices ) );
    std::stringstream msoftBuilder;
    msoftBuilder << "MSOFT[" << msoftIndex << "]";
    size_t const linearMassIndex( RegisterBlockEntry( msoftBuilder.str() ) );
    return AddNewDerivedParameter( caseString,
                                   new SlhaMassSquaredDiagonalFunctionoid(
                                              numberOfDistinctActiveParameters,
                                                               squareMassIndex,
                                                           linearMassIndex ) );
  }


  // This creates a temporary functionoid and passes it parameter values from
  // OnceOff(...) to return the appropriate value for caseString.
  inline double
  SlhaBlocksWithSpecialCasesManager::OnceOffSlhaOneOrTwoCompatibleMassSquared(
                                                 std::string const& caseString,
                                                       char const sfermionType,
                                                    char const generationIndex,
                                                 unsigned int const msoftIndex,
                                          double const logarithmOfScale ) const
  {
    std::string blockNameWithIndices( "MS?2[?,?]" );
    blockNameWithIndices[ 2 ] = sfermionType;
    blockNameWithIndices[ 5 ] = blockNameWithIndices[ 7 ] = generationIndex;
    double const squareMass( OnceOffParameter( blockNameWithIndices,
                                               logarithmOfScale ) );
    std::stringstream msoftBuilder;
    msoftBuilder << "MSOFT[" << msoftIndex << "]";
    double const linearMass( OnceOffParameter( msoftBuilder.str(),
                                               logarithmOfScale ) );
    SlhaMassSquaredDiagonalFunctionoid temporaryParameter( 0,
                                                           0,
                                                           0 );
    return temporaryParameter( squareMass,
                               linearMass );
  }

  // This returns a string which is the concatenated set of strings from
  // parameter functionoids giving their Python evaluations.
  inline std::string
  SlhaBlocksWithSpecialCasesManager::ParametersInPythonFunction(
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

} /* namespace VevaciousPlusPlus */

#endif /* SLHABLOCKSWITHSPECIALCASESMANAGER_HPP_ */
