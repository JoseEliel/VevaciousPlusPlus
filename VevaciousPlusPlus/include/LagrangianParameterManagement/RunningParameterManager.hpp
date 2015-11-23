/*
 * RunningParameterManager.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef RUNNINGPARAMETERMANAGER_HPP_
#define RUNNINGPARAMETERMANAGER_HPP_

#include "../LagrangianParameterManagement/SlhaManager.hpp"
#include "../LagrangianParameterManagement/ParameterUpdatePropagator.hpp"
#include "CommonIncludes.hpp"
#include "../PotentialEvaluation/ParameterFunctionoid.hpp"
#include "../PotentialEvaluation/ParameterFunctionoids/BinaryOperationFunctionoid.hpp"
#include "../PotentialEvaluation/ParameterFunctionoids/ConstantFunctionoid.hpp"
#include "../PotentialEvaluation/ParameterFunctionoids/SlhaFunctionoid.hpp"
#include "../PotentialEvaluation/ParameterFunctionoids/UnaryOperationFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class RunningParameterManager : public SlhaManager,
                                  public ParameterUpdatePropagator
  {
  public:
    typedef LHPC::SLHA::SparseManyIndexedBlock< double > RestrictedSlhaBlock;

    // This puts all variables with index brackets into a consistent form.
    virtual std::string
    FormatVariable( std::string const& unformattedVariable );


    RunningParameterManager(
                      LagrangianParameterManager& lagrangianParameterManager );
    virtual ~RunningParameterManager();


    // This adds to the list of SLHA blocks that will be read in when updating
    // running parameters. An alias for the block can be provided.
    void AddValidSlhaBlock( std::string slhaBlock,
                            std::string const& aliasString );

    // This creates a new ParameterFunctionoid and a string mapping to it.
    void CreateDerivedParameter( std::string const& parameterString );

    // This returns a pointer to the ParameterFunctionoid mapped to by
    // parameterString, returning NULL if there is none.
    ParameterFunctionoid* GetFunctionoid( std::string const& parameterString );

    // This updates all SlhaFunctionoids with data from a new SLHA file, sets
    // lowestBlockScale and highestBlockScale, then updates the functionoids in
    // parameterFunctionoidPointers to take their values at lowestBlockScale.
    virtual void ReadFile( std::string const& slhaFilename );

    // This returns the lowest Q found among the SLHA blocks. It returns 0.0
    // if no block had a scale greater than 0.0 GeV.
    double LowestBlockScale() const{ return lowestBlockScale; }

    // This returns the highest Q found among the SLHA blocks. It returns 0.0
    // if no block had a scale greater than 0.0 GeV.
    double HighestBlockScale() const{ return highestBlockScale; }

    // This returns a string that should be valid Python to set the running
    // parameters in order for a given renormalization scale.
    std::string RunningParametersAsPython() const;

    // This is just a dummy to make things work while I have it side-by-side to
    // the new version for comparison.
    virtual void UpdateSelfForNewParameterPoint(
                 LagrangianParameterManager const& lagrangianParameterManager )
    { ; /* This does nothing. */ }

    // This is just a dummy to make things work while I have it side-by-side to
    // the new version for comparison.
    virtual double OnceOffParameter( std::string const& parameterName,
                                     double const logarithmOfScale ) const
    { return 0.0; }

    // This is just a dummy to make things work while I have it side-by-side to
    // the new version for comparison.
    virtual std::pair< bool, size_t >
    RegisterParameter( std::string const& parameterName )
    { return std::make_pair( false, -1 ); }

    // This is just a dummy to make things work while I have it side-by-side to
    // the new version for comparison.
    virtual void ParameterValues( double logarithmOfScale,
                               std::vector< double >& destinationVector ) const
    { ; /* This does nothing. */ }

    // This returns a scale which is appropriate for using for a fixed-scale
    // calculation for the current parameter point.
    virtual double MinimumEvaluationScale() const
    { return lowestBlockScale; }

    // This returns a scale which is appropriate for using for a fixed-scale
    // calculation for the current parameter point.
    virtual double AppropriateSingleFixedScale() const
    { return lowestBlockScale; }

    // This returns a scale which is appropriate for using for a fixed-scale
    // calculation for the current parameter point.
    virtual double MaximumEvaluationScale() const
    { return highestBlockScale; }


  protected:
    // This puts all index brackets into a consistent form.
    virtual std::string
    FormatIndexBracketContent( std::string const& unformattedBracketContent );

    std::vector< ParameterFunctionoid* > parameterFunctionoidPointers;
    std::vector< SlhaFunctionoid* > slhaFunctionoidPointers;
    std::map< std::string, ParameterFunctionoid* > parameterFunctionoidMap;
    std::map< std::string, ParameterFunctionoid* >::iterator parameterFinder;
    std::map< double, ParameterFunctionoid* > constantFunctionoidMap;
    std::map< double, ParameterFunctionoid* >::iterator constantFinder;
    LHPC::SlhaParser slhaParser;
    std::vector< RestrictedSlhaBlock* > slhaBlockPointers;
    std::map< std::string, RestrictedSlhaBlock* > slhaBlockMap;
    std::map< std::string, RestrictedSlhaBlock* >::iterator slhaBlockFinder;
    std::map< std::string, std::string > slhaAliasMap;
    std::map< std::string, std::string >::iterator slhaAliasFinder;
    double lowestBlockScale;
    double highestBlockScale;

    // This returns a pointer to the ParameterFunctionoid for the block named
    // blockName with indices given by indexString, making a new one if there
    // is none already. If the block name was invalid though, nothing is made
    // and NULL is returned.
    ParameterFunctionoid* GetSlhaFunctionoid( std::string const& blockName,
                                              std::string const& indexString,
                                          std::string const& parameterString );

    // This returns a pair of functionoid pointers which match the arguments
    // in commaSeparatedAliases.
    std::pair< ParameterFunctionoid*, ParameterFunctionoid* >
    FindFunctionoidPair( std::string const& commaSeparatedAliases );

    // This makes a name based on how many functionoids are already in
    // parameterFunctionoidPointers.
    std::string FunctionoidName() const;

    // This is just a dummy to make things work while I have it side-by-side to
    // the new version for comparison.
    virtual void PrepareNewParameterPoint( std::string const& newInput )
    { ; /* This does nothing. */ }
  };




  // This puts all variables with index brackets into a consistent form.
  inline std::string RunningParameterManager::FormatVariable(
                                       std::string const& unformattedVariable )
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
  inline std::string RunningParameterManager::FormatIndexBracketContent(
                                 std::string const& unformattedBracketContent )
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

  // This makes a name based on how many functionoids are already in
  // parameterFunctionoidPointers.
  inline std::string RunningParameterManager::FunctionoidName() const
  {
    std::stringstream stringBuilder;
    stringBuilder << "rp" << parameterFunctionoidPointers.size();
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* RUNNINGPARAMETERMANAGER_HPP_ */
