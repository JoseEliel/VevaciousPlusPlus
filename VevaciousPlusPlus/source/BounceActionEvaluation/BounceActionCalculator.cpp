/*
 * BounceActionCalculator.cpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionCalculator.hpp"

namespace VevaciousPlusPlus
{

  BounceActionCalculator::BounceActionCalculator()
  {
    // This does nothing.
  }

  BounceActionCalculator::~BounceActionCalculator()
  {
    // This does nothing.
  }


  // This plots the fields as functions of the spatial variables in a file
  // called plotFilename in .eps format, with each field plotted in the color
  // given by fieldColors: the field with index i is plotted in the color
  // given by fieldColors[ i ]. An empty string indicates that the field
  // should not be plotted.
  void BounceActionCalculator::PlotBounceConfiguration(
                                                  TunnelPath const& tunnelPath,
                                            BubbleProfile const& bubbleProfile,
                                  std::vector< std::string > const& fieldNames,
                                 std::vector< std::string > const& fieldColors,
                                               std::string const& plotFilename,
                                     unsigned int const plotResolution ) const
  {
    BOL::TwoDimensionalDataPlotter bubblePlotter( "/usr/bin/gnuplot",
                                                  plotFilename );
    double const radialStepSize( bubbleProfile.MaximumPlotRadius()
                                 / static_cast< double >( plotResolution ) );
    size_t numberOfPlottedFields( std::min( fieldColors.size(),
                                            tunnelPath.NumberOfFields() ) );
    BOL::TwoDimensionalDataPlotter::PlotDataVector plotData;
    BOL::TwoDimensionalDataPlotter::DoublePairVectorWithStringPair fieldData;
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfPlottedFields;
         ++fieldIndex )
    {
      if( !(fieldColors[ fieldIndex ].empty()) )
      {
        fieldData.second.first.assign( fieldColors[ fieldIndex ] );
        fieldData.second.second.assign( fieldNames[ fieldIndex ] );
        plotData.push_back( fieldData );
      }
    }
    numberOfPlottedFields = plotData.size();
    std::vector< double >
    fieldConfiguration( tunnelPath.NumberOfFields() );
    std::cout << "Bubble profile:" << std::endl;
    for( unsigned int plotIndex( 0 );
         plotIndex < plotResolution;
         ++plotIndex )
    {
      double const radialValue( plotIndex * radialStepSize );
      double const auxiliaryValue( bubbleProfile.AuxiliaryAt( radialValue ) );
      double const
      auxiliarySlope( bubbleProfile.AuxiliarySlopeAt( radialValue ) );
      std::cout << "r = " << radialValue << ", p = " << auxiliaryValue
      << ", dp/dr = " << auxiliarySlope << ", "
      << tunnelPath.FieldsString( auxiliaryValue )
      << std::endl;
      tunnelPath.PutOnPathAt( fieldConfiguration,
                              auxiliaryValue );
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfPlottedFields;
           ++fieldIndex )
      {
        plotData[ fieldIndex ].first.push_back( std::make_pair( radialValue,
                                          fieldConfiguration[ fieldIndex ] ) );
      }
    }
    std::cout << std::endl;
    bubblePlotter.plotData( plotData,
                            "rho/(1/GeV)",
                            "field/GeV" );
  }

} /* namespace VevaciousPlusPlus */
