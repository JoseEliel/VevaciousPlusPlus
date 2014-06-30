/*
 * PotentialMinimum.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POTENTIALMINIMUM_HPP_
#define POTENTIALMINIMUM_HPP_

#include "CommonIncludes.hpp"
#include "MinuitMinimization/MinuitMinimum.hpp"

namespace VevaciousPlusPlus
{
  class PotentialMinimum : public MinuitMinimum
  {
  public:
    PotentialMinimum( std::vector< double > const& fieldConfiguration,
                      double const potentialDepth );
    PotentialMinimum( MinuitMinimum const& minuitMinimum );
    PotentialMinimum();
    virtual
    ~PotentialMinimum();


    // This returns the sum of the squares of the differences in the field
    // values of this PotentialMinimum with comparisonFields.
    double
    SquareDistanceTo( std::vector< double > const& comparisonFields ) const;

    // This returns the sum of the squares of the differences in the field
    // values of this PotentialMinimum with comparisonMinimum.
    double SquareDistanceTo( PotentialMinimum const& comparisonMinimum ) const
    { return SquareDistanceTo( comparisonMinimum.FieldConfiguration() ); }

    // This returns the sum of the squares of the field values.
    double LengthSquared() const;

    std::vector< double > const& FieldConfiguration() const
    { return variableValues; }

    double PotentialValue() const{ return functionValue; }


    // This prints the minimum as an empty XML element.
    std::string AsEmptyXmlElement( std::string const& elementName,
                          std::vector< std::string > const& fieldNames ) const;

    // This prints the minimum as an XML element in the style chosen for
    // Vevacious output.
    std::string AsVevaciousXmlElement( std::string const& elementName,
                          std::vector< std::string > const& fieldNames ) const;

    // This prints the minimum in a form that Mathematica can understand.
    std::string
    AsMathematica( std::vector< std::string > const& fieldNames ) const;
  };




  // This returns the sum of the squares of the differences in the field
  // values of this PotentialMinimum with comparisonFields.
  inline double PotentialMinimum::SquareDistanceTo(
                          std::vector< double > const& comparisonFields ) const
  {
    double returnDouble( 0.0 );
    double fieldDifference( 0.0 );
    for( unsigned int fieldIndex( 0 );
         fieldIndex < variableValues.size();
         ++fieldIndex )
    {
      fieldDifference
      = ( variableValues[ fieldIndex ] - comparisonFields[ fieldIndex ] );
      returnDouble += ( fieldDifference * fieldDifference );
    }
    return returnDouble;
  }

  // This returns the sum of the squares of the field values.
  inline double PotentialMinimum::LengthSquared() const
  {
    double returnDouble( 0.0 );
    for( unsigned int fieldIndex( 0 );
         fieldIndex < variableValues.size();
         ++fieldIndex )
    {
      returnDouble += ( variableValues[ fieldIndex ]
                       * variableValues[ fieldIndex ] );
    }
    return returnDouble;
  }


  // This prints the minimum as an empty XML element.
  inline std::string
  PotentialMinimum::AsEmptyXmlElement( std::string const& elementName,
                           std::vector< std::string > const& fieldNames ) const
  {
    std::stringstream stringBuilder;
    stringBuilder << "<" << elementName;
    for( unsigned int fieldIndex( 0 );
         fieldIndex < variableValues.size();
         ++fieldIndex )
    {
      stringBuilder << " " << fieldNames[ fieldIndex ] << "=\""
      << variableValues[ fieldIndex ] << "\"";
    }
    stringBuilder << " PotentialDepth=\"" << functionValue << "\" />";
    return stringBuilder.str();
  }

  // This prints the minimum as an XML element in the style chosen for
  // Vevacious output.
  inline std::string
  PotentialMinimum::AsVevaciousXmlElement( std::string const& elementName,
                           std::vector< std::string > const& fieldNames ) const
  {
    std::stringstream stringBuilder;
    stringBuilder << "  <" << elementName << ">\n"
    "    <FieldValues>\n";
    for( unsigned int fieldIndex( 0 );
         fieldIndex < variableValues.size();
         ++fieldIndex )
    {
      stringBuilder << "      " << variableValues[ fieldIndex ] << " <!-- "
      << fieldNames[ fieldIndex ] << " -->\n";
    }
    stringBuilder << "    </FieldValues>\n"
    "    <RelativeDepth>\n"
    "      " << functionValue
    << " <!-- Potential in GeV^4 minus value for all fields = 0 -->\n"
    "    </RelativeDepth>\n"
    "  </" << elementName << ">\n";
    return stringBuilder.str();
  }

  // This prints the minimum in a form that Mathematica can understand.
  inline std::string PotentialMinimum::AsMathematica(
                           std::vector< std::string > const& fieldNames ) const
  {
    std::stringstream stringBuilder;
    stringBuilder << "{ {";
    for( unsigned int fieldIndex( 0 );
         fieldIndex < variableValues.size();
         ++fieldIndex )
    {
      if( fieldIndex != 0 )
      {
        stringBuilder << ",";
      }
      stringBuilder << " " << fieldNames[ fieldIndex ] << " -> "
      << variableValues[ fieldIndex ];
    }
    stringBuilder << " }, PotentialDepth -> " << functionValue << " }";
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* POTENTIALMINIMUM_HPP_ */
