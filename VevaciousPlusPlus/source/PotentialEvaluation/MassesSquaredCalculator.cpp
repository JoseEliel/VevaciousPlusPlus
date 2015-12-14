/*
 * MassesSquaredCalculator.cpp
 *
 *  Created on: Oct 09, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialEvaluation/MassesSquaredCalculator.hpp"

namespace VevaciousPlusPlus
{

  MassesSquaredCalculator::MassesSquaredCalculator(
                   std::map< std::string, std::string > const& attributeMap ) :
    multiplicityFactor( 1.0 ),
    spinType( notSet )
  {
    std::map< std::string, std::string >::const_iterator
    attributeFinder( attributeMap.find( "MultiplicityFactor" ) );
    if( attributeFinder != attributeMap.end() )
    {
      multiplicityFactor
      = LHPC::ParsingUtilities::StringToDouble( attributeFinder->second );
    }
    attributeFinder = attributeMap.find( "SpinType" );
    if( attributeFinder != attributeMap.end() )
    {
      if( attributeFinder->second == "ScalarBoson" )
      {
        spinType = scalarBoson;
      }
      else if( attributeFinder->second == "WeylFermion" )
      {
        spinType = weylFermion;
      }
      else if( attributeFinder->second == "GaugeBoson" )
      {
        spinType = gaugeBoson;
      }
    }
  }

  MassesSquaredCalculator::MassesSquaredCalculator(
                                  MassesSquaredCalculator const& copySource ) :
    multiplicityFactor( copySource.multiplicityFactor ),
    spinType( copySource.spinType )
  {
    // This constructor is just an initialization list.
  }

  MassesSquaredCalculator::MassesSquaredCalculator() :
    multiplicityFactor( 0.0 ),
    spinType( notSet )
  {
    // This constructor is just an initialization list.
  }

  MassesSquaredCalculator::~MassesSquaredCalculator()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
