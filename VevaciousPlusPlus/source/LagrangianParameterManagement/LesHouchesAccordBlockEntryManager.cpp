/*
 * LesHouchesAccordBlockEntryManager.cpp
 *
 *  Created on: Oct 29, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/LesHouchesAccordBlockEntryManager.hpp"

namespace VevaciousPlusPlus
{

  LesHouchesAccordBlockEntryManager::LesHouchesAccordBlockEntryManager(
                                 std::set< std::string > const& validBlocks ) :
    LagrangianParameterManager(),
    numberOfDistinctActiveParameters( 0 ),
    activeParametersToIndices(),
    activeInterpolatedParameters(),
    validBlocks( validBlocks ),
    lhaParser()
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "LesHouchesAccordBlockEntryManager using SlhaPolynomialFitBlockEntry"
        " objects! Probably should change this to"
        " SlhaLinearlyInterpolatedBlockEntry objects.";
    std::cout << std::endl;/**/
  }

  LesHouchesAccordBlockEntryManager::~LesHouchesAccordBlockEntryManager()
  {
    // This does nothing.
  }


  // This checks to see if the parameter name has already been registered,
  // and if so, returns true paired with the index to its functionoid. If
  // not, then it checks to see if the parameter name corresponds to an entry
  // in a valid LHA block, and if so, adds an interpolation functionoid to
  // the set of functionoids and returns true paired with the index of the
  // new functionoid. If the parameter name does not fall into any of the
  // above cases, then false is returned, paired with -1, rolling over to the
  // maximum value of size_t as it is unsigned.
  std::pair< bool, size_t >
  LesHouchesAccordBlockEntryManager::RegisterParameter(
                                             std::string const& parameterName )
  {
    std::map< std::string, size_t >::const_iterator
    alreadyExistsResult( activeParametersToIndices.find( parameterName ) );
    if( alreadyExistsResult != activeParametersToIndices.end() )
    {
      return std::pair< bool, size_t >( true, alreadyExistsResult->second );
    }

    // If it has not already been registered, we check whether it is a valid
    // LHA block parameter which has not yet been registered.
    if( RefersToValidBlock( parameterName ) )
    {
      activeInterpolatedParameters.push_back( LhaBlockEntryInterpolator(
                                              numberOfDistinctActiveParameters,
                                                                     lhaParser,
                                                             parameterName ) );
      // We tersely update numberOfDistinctActiveParameters with
      // post-increment ++ so that the pair has the correct index.
      return std::pair< bool, size_t >( true,
                                        numberOfDistinctActiveParameters++ );
    }

    // If parameterName did not correspond to an existing parameter or a
    // valid new LHA block entry parameter, then we return false with the
    // maximum size_t, which should cause out-of-range exceptions if used
    // despite the false.
    return std::pair< bool, size_t >( false,
                                      -1 );
  }

} /* namespace VevaciousPlusPlus */
