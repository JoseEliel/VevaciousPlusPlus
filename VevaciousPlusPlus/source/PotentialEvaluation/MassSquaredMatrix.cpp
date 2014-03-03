/*
 * MassSquaredMatrix.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  MassSquaredMatrix::MassSquaredMatrix(
                   std::map< std::string, std::string > const& attributeMap ) :
    matrixElements(),
    eigenMatrix(),
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
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MassSquaredMatrix::MassSquaredMatrix( ... )";
    std::cout << std::endl;/**/
  }

  MassSquaredMatrix::~MassSquaredMatrix()
  {
    // This does nothing.
  }


  // This returns the eigenvalues of the matrix.
  std::vector< double > const& MassSquaredMatrix::MassesSquared(
                              std::vector< double > const& fieldConfiguration )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MassSquaredMatrix::MassesSquared( ... )";
    std::cout << std::endl;

    return massesSquared;/**/
  }

} /* namespace VevaciousPlusPlus */
