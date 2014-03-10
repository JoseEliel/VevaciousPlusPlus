/*
 * MassesSquaredFromPolynomials.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  MassesSquaredFromPolynomials::MassesSquaredFromPolynomials(
                   std::map< std::string, std::string > const& attributeMap ) :
    massesSquared(),
    multiplicityFactor( 1.0 ),
    spinType( notSet )
  {
    std::map< std::string, std::string >::const_iterator
    attributeFinder( attributeMap.find( "MultiplicityFactor" ) );
    if( attributeFinder != attributeMap.end() )
    {
      multiplicityFactor
      = BOL::StringParser::stringToDouble( attributeFinder->second );
    }
    attributeFinder = attributeMap.find( "SpinType" );
    if( attributeFinder != attributeMap.end() )
    {
      if( attributeFinder->second.compare( "ScalarBoson" ) == 0 )
      {
        spinType = scalarBoson;
      }
      else if( attributeFinder->second.compare( "WeylFermion" ) == 0 )
      {
        spinType = weylFermion;
      }
      else if( attributeFinder->second.compare( "GaugeBoson" ) == 0 )
      {
        spinType = gaugeBoson;
      }
    }
  }

  MassesSquaredFromPolynomials::MassesSquaredFromPolynomials(
                             MassesSquaredFromPolynomials const& copySource ) :
    massesSquared( copySource.massesSquared ),
    multiplicityFactor( copySource.multiplicityFactor ),
    spinType( copySource.spinType )
  {
    // This constructor is just an initialization list.
  }

  MassesSquaredFromPolynomials::MassesSquaredFromPolynomials() :
    massesSquared(),
    multiplicityFactor( NAN ),
    spinType( notSet )
  {
    // This constructor is just an initialization list.
  }

  MassesSquaredFromPolynomials::~MassesSquaredFromPolynomials()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
