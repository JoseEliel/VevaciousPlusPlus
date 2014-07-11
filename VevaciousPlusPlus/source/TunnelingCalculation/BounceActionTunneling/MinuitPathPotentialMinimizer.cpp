/*
 * MinuitPathPotentialMinimizer.cpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "TunnelingCalculation/BounceActionTunneling/MinuitPathPotentialMinimizer.hpp"

namespace VevaciousPlusPlus
{

  MinuitPathPotentialMinimizer::MinuitPathPotentialMinimizer()
  {
    BOL::AsciiXmlParser outerArgumentParser;
    BOL::AsciiXmlParser innerArgumentParser;
    outerArgumentParser.loadString( xmlArguments );
    std::string pathFactoryType( "PathFromNodesOnParallelPlanes" );
    std::string pathFactoryArguments( "" );
    while( outerArgumentParser.readNextElement() )
    {
      if( outerArgumentParser.currentElementNameMatches(
                                                     "PathParameterization" ) )
      {
        innerArgumentParser.loadString(
                       outerArgumentParser.getTrimmedCurrentElementContent() );
        innerArgumentParser.readNextElement();
        pathFactoryType.assign(
                       innerArgumentParser.getTrimmedCurrentElementContent() );
        innerArgumentParser.readNextElement();
        pathFactoryArguments.assign(
                       innerArgumentParser.getTrimmedCurrentElementContent() );
      }
      else if( outerArgumentParser.currentElementNameMatches(
                                            "NumberOfPotentialSamplePoints" ) )
      {
        minuitToleranceFraction
        = BOL::StringParser::stringToDouble(
                       outerArgumentParser.getTrimmedCurrentElementContent() );
      }
      else if( outerArgumentParser.currentElementNameMatches(
                                                           "MinuitStrategy" ) )
      {
        minuitStrategy
        = BOL::StringParser::stringToInt(
                       outerArgumentParser.getTrimmedCurrentElementContent() );
      }
      else if( outerArgumentParser.currentElementNameMatches(
                                                          "MinuitTolerance" ) )
      {
        minuitToleranceFraction
        = BOL::StringParser::stringToDouble(
                       outerArgumentParser.getTrimmedCurrentElementContent() );
      }
    }
    if( pathFactoryType.compare( "LinearSplineThroughNodes" ) == 0 )
    {
      pathFactory
      = new LinearSplineThroughNodesFactory( pathFactoryArguments );
    }
    else if( pathFactoryType.compare( "QuadraticSplineThroughNodes" ) == 0 )
    {
      pathFactory
      = new QuadraticSplineThroughNodesFactory( pathFactoryArguments );
    }
    else if( pathFactoryType.compare( "PolynomialThroughNodes" ) == 0 )
    {
      pathFactory = new PolynomialThroughNodesFactory( pathFactoryArguments );
    }
    pathNodes = &(pathFactory->NodesFromParameterization());
  }

  MinuitPathPotentialMinimizer::~MinuitPathPotentialMinimizer()
  {
    // TODO Auto-generated destructor stub
  }

} /* namespace VevaciousPlusPlus */
