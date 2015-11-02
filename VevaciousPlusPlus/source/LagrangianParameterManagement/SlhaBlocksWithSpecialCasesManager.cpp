/*
 * SlhaBlocksWithSpecialCasesManager.cpp
 *
 *  Created on: Nov 2, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/SlhaBlocksWithSpecialCasesManager.hpp"

namespace VevaciousPlusPlus
{

  SlhaBlocksWithSpecialCasesManager::SlhaBlocksWithSpecialCasesManager()
  {
    // TODO Auto-generated constructor stub

  }

  SlhaBlocksWithSpecialCasesManager::~SlhaBlocksWithSpecialCasesManager()
  {
    for( size_t deletionIndex( 0 );
         deletionIndex < activeDerivedParameters.size();
         ++deletionIndex )
    {
      delete activeDerivedParameters[ deletionIndex ];
    }
  }

} /* namespace VevaciousPlusPlus */
