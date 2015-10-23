/*
 * SlhaInterpolatedParameterManager.cpp
 *
 *  Created on: Oct 23, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/SlhaInterpolatedParameterManager.hpp"

namespace VevaciousPlusPlus
{

  SlhaInterpolatedParameterManager::SlhaInterpolatedParameterManager()
  {
    // TODO Auto-generated constructor stub
  }

  SlhaInterpolatedParameterManager::~SlhaInterpolatedParameterManager()
  {
    // TODO Auto-generated destructor stub
  }


  // This checks to see if the parameter name has already been registered,
  // and if so, returns true paired with the index to its functionoid. If
  // not, it checks to see if it corresponds to a valid special case, and if
  // so, it adds the special case functionoid to the set of active
  // functionoids and returns true paired with the index of the newly active
  // functionoid. If not, then it checks to see if the parameter name
  // corresponds to an entry in a valid SLHA block, and if so, adds the
  // normal SLHA functionoid to the set of active functionoids and returns
  // true paired with the index of the newly active functionoid. If the
  // parameter name does not fall into any of the above cases, then false is
  // returned, paired with -1, rolling over to the maximum value of size_t as
  // it is unsigned.
  std::pair< bool, size_t >
  SlhaInterpolatedParameterManager::RegisterParameter(
                                             std::string const& parameterName )
  {
    std::map< std::string, size_t >::const_iterator
    alreadyExistsResult( activeParametersToIndices.find( parameterName ) );
    if( alreadyExistsResult != activeParametersToIndices.end() )
    {
      return std::pair< bool, size_t >( true, alreadyExistsResult->second );
    }

    // Next we have to find if the parameter name is a special case.
    std::map< std::string, SlhaDerivedFunctionoid* >::const_iterator
    specialCaseResult( validSpecialCases.find( parameterName ) );
    if( specialCaseResult != validSpecialCases.end() )
    {
      SlhaDerivedFunctionoid* foundFunctionoid = specialCaseResult->second;

      // No matter if the functionoid already exists with an active index
      // through having already been registered with a different alias
      // (several names might map to the same derived functionoid, such as
      // TE[2,2] and Te22), this name now maps to the found functionoid.
      activeParametersToIndices[ parameterName ] = foundFunctionoid;
      if( !(foundFunctionoid->IsActive()) )
      {
        // If the functionoid has not already been registered, it is added to
        // the list of those which will be updated by
        // PrepareNewParameterPoint, and given the next appropriate index
        // equal to the number of existing active parameters, and this index
        // is incremented afterwards.
        activeSpecialCases.push_back( foundFunctionoid );
        foundFunctionoid->MakeActive( numberOfDistinctActiveParameters++ );
      }
      return std::pair< bool, size_t >( true,
                                   foundFunctionoid->IndexInValuesVector() );
    }

    // Finally we check to find if it is a valid normal SLHA parameter which
    // has not yet been registered.
    std::string indexString;
    LHPC::SLHA::SparseManyIndexedBlock< double > const*
    slhaBlock( BlockFromParameterName( parameterName,
                                       indexString ) );
    if( slhaBlock != NULL )
    {
      activeNormalSlhaParameters.push_back(
                                        SlhaInterpolatedParameterFunctionoid(
                                            numberOfDistinctActiveParameters,
                                                                   slhaBlock,
                                                             indexString ) );
      // We tersely update numberOfDistinctActiveParameters with
      // post-increment ++ so that the pair has the correct index.
      return std::pair< bool, size_t >( true,
                                        numberOfDistinctActiveParameters++ );
    }

    // If parameterName did not correspond to an existing parameter, a valid
    // new special case, or a valid new normal SLHA parameter, then we return
    // false with the maximum size_t, which should cause out-of-range
    // exceptions if used despite the false.
    return std::pair< bool, size_t >( false,
                                      -1 );
  }

} /* namespace VevaciousPlusPlus */
