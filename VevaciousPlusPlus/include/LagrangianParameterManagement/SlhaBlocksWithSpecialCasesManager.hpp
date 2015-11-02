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
#include "SlhaSourcedParameterFunctionoid.hpp"
#include "SlhaDerivedFunctionoids/SlhaDsbHiggsVevFunctionoid.hpp"
#include "SlhaDerivedFunctionoids/SlhaHiggsMixingBilinearFunctionoid.hpp"
#include "SlhaDerivedFunctionoids/SlhaMassSquaredDiagonalFunctionoid.hpp"
#include "SlhaDerivedFunctionoids/SlhaTrilinearDiagonalFunctionoid.hpp"
#include "SlhaDerivedFunctionoids/SlhaTwoSourceFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class SlhaBlocksWithSpecialCasesManager :
                                       public LesHouchesAccordBlockEntryManager
  {
  public:
    SlhaBlocksWithSpecialCasesManager( std::string const& validBlocksString );
    SlhaBlocksWithSpecialCasesManager(
                                   std::set< std::string > const& validBlocks);
    virtual ~SlhaBlocksWithSpecialCasesManager();


  protected:
    std::vector< SlhaSourcedParameterFunctionoid* > activeDerivedParameters;


    // This adds newParameter to activeDerivedParameters and updates
    // numberOfDistinctActiveParameters and activeParametersToIndices, and
    // returns true paired with the index of newParameter. (The index mapped to
    // in activeParametersToIndices and used as part of the return value is
    // based on numberOfDistinctActiveParameters before it gets incremented,
    // rather than using newParameter->IndexInValuesVector().)
    std::pair< bool, size_t >
    AddNewDerivedParameter( std::string const& parameterName,
                            SlhaSourcedParameterFunctionoid* newParameter );

    // This checks parameterName against all the special cases.
    virtual std::pair< bool, size_t >
    RegisterUnregisteredParameter( std::string const& parameterName )
    {
      // If the parameter does not already exist, we create its functionoid and
      // add all its aliases to aliasesToActiveParameters, add the functionoid
      // to activeDerivedParameters, and update
      switch( parameterName )
      {
        /*
        vd
        vu
        mu + tree
        Bmu + tree
        mHdSq +  tree
        mHuSq +  tree
        Te11, Td11, Tu11, etc.
        M2L11, M2E11, M2Q11, M2U11, M2D11, etc.
        */


        case "DsbVd":
        case "DsbVu":
          SlhaSourcedParameterFunctionoid const&
          vevLength( RegisterBlockEntry( "HMIX[ 3 ]" ) );
          SlhaSourcedParameterFunctionoid const&
          tanBeta( RegisterBlockEntry( "HMIX[ 2 ]" ) );

          return AddNewDerivedParameter( parameterName,
              new SlhaDsbHiggsVevFunctionoid( numberOfDistinctActiveParameters,
                                              vevLength,
                                              tanBeta,
                                  ( parameterName.compare( "DsbVu" ) == 0 ) ));
          break;
        case "Bmu":
        case "m3Sq":
          SlhaSourcedParameterFunctionoid const&
          treePseudoscalarMassSquared( RegisterBlockEntry( "HMIX[ 4 ]" ) );
          SlhaSourcedParameterFunctionoid const&
          tanBeta( RegisterBlockEntry( "HMIX[ 2 ]" ) );
          // We ensure that both aliases map to this parameter while
          // numberOfDistinctActiveParameters is still equal to the index it
          // will have.

          activeParametersToIndices[ "Bmu" ]
          = activeParametersToIndices[ "m3Sq" ]
          = numberOfDistinctActiveParameters;
          return AddNewDerivedParameter( parameterName,
                                        new SlhaHiggsMixingBilinearFunctionoid(
                                              numberOfDistinctActiveParameters,
                                                   treePseudoscalarMassSquared,
                                                                   tanBeta ) );
        default:
          return std::pair< bool, size_t >( false,
                                            -1 );
      }
    }
  };

  // This adds newParameter to activeDerivedParameters and updates
  // numberOfDistinctActiveParameters and activeParametersToIndices, and
  // returns true paired with the index of newParameter. (The index mapped to
  // in activeParametersToIndices and used as part of the return value is based
  // on numberOfDistinctActiveParameters before it gets incremented, rather
  // than using newParameter->IndexInValuesVector().)
  inline std::pair< bool, size_t >
  SlhaBlocksWithSpecialCasesManager::AddNewDerivedParameter(
                                              std::string const& parameterName,
                                SlhaSourcedParameterFunctionoid* newParameter )
  {
    activeDerivedParameters.push_back( newParameter );
    activeParametersToIndices[ parameterName ]
    = numberOfDistinctActiveParameters;
    // We tersely update numberOfDistinctActiveParameters with
    // post-increment ++ so that the pair has the correct index.
    return std::pair< bool, size_t >( true,
                                      numberOfDistinctActiveParameters++ );
  }

} /* namespace VevaciousPlusPlus */

#endif /* SLHABLOCKSWITHSPECIALCASESMANAGER_HPP_ */
