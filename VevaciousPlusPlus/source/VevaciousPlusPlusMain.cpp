/*
 * VevaciousPlusPlusMain.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 */

#include "VevaciousPlusPlus.hpp"
#include "LHPC/Utilities/RestrictedXmlParser.hpp"
#include "Utilities/FilePlaceholderManager.hpp"


int main( int argumentCount,
          char** argumentCharArrays )
{
  if( argumentCount != 2 )
  {
    std::cout
    << std::endl
    << "VevaciousPlusPlus requires a single argument: the name of the XML"
    << " file which directs the code to further initialization files.";
    std::cout << std::endl;
    std::cout << "An example \"main input\" file is provided:"
    << " \"bin/VevaciousPlusPlusMainInput.xml\", which has comments describing"
    << " each of the elements.";
    std::cout << std::endl;
  }
  else
  {
    std::string inputFilename( argumentCharArrays[ 1 ] );
    std::string initializationFile( "" );
    std::vector< std::pair< std::string, std::string > > parameterPoints;
    LHPC::RestrictedXmlParser xmlParser;
    xmlParser.OpenRootElementOfFile( inputFilename );
    while( xmlParser.ReadNextElement() )
    {
      if( xmlParser.CurrentName() == "InitializationFile" )
      {
        initializationFile = xmlParser.TrimmedCurrentBody();
      }
      else if( ( xmlParser.CurrentName() == "SingleParameterPoint" )
               ||
               ( xmlParser.CurrentName() == "ParameterPointSet" ) )
      {
        parameterPoints.push_back( std::make_pair( xmlParser.CurrentName(),
                                                   xmlParser.CurrentBody() ) );
      }
    }
    xmlParser.CloseFile();

    // The default is to construct the VevaciousPlusPlus object with an input
    // initialization file name, which will create a PotentialMinimizer and
    // TunnelingCalculator internal to the VevaciousPlusPlus object, managing
    // the memory allocated in doing so, with all components initialized
    // according to the XML input file. Alternatively, one can create the
    // PotentialMinimizer and TunnelingCalculator components externally and
    // pass them to the other constructor.
    VevaciousPlusPlus::VevaciousPlusPlus
    vevaciousPlusPlus( initializationFile );

    std::string runPointInput( "" );
    std::string outputFilename( "" );
    bool appendLhaOutputToLhaInput( false );
    std::string inputFolder( "" );
    std::string outputFolder( "" );
    for( std::vector< std::pair< std::string, std::string > >::const_iterator
         parameterElement( parameterPoints.begin() );
         parameterElement != parameterPoints.end();
         ++parameterElement )
    {
      appendLhaOutputToLhaInput = false;
      xmlParser.LoadString( parameterElement->second );
      if( parameterElement->first == "SingleParameterPoint" )
      {
        while( xmlParser.ReadNextElement() )
        {
          if( xmlParser.CurrentName() == "RunPointInput" )
          {
            runPointInput = xmlParser.TrimmedCurrentBody();
          }
          else if( xmlParser.CurrentName() == "OutputFilename" )
          {
            outputFilename = xmlParser.TrimmedCurrentBody();
          }
          else if( xmlParser.CurrentName() == "AppendLhaOutputToLhaInput" )
          {
            appendLhaOutputToLhaInput = true;
          }
        }

        if( !appendLhaOutputToLhaInput && outputFilename.empty() )
        {
          std::stringstream errorBuilder;
          errorBuilder << "Single parameter point cannot be run without either"
          << " an output XML file name or the <AppendLhaOutputToLhaInput />"
          << " option.";
          throw std::runtime_error( errorBuilder.str() );
        }
        vevaciousPlusPlus.RunPoint( runPointInput );
        if( !(outputFilename.empty()) )
        {
          vevaciousPlusPlus.WriteResultsAsXmlFile( outputFilename );
        }
        if( appendLhaOutputToLhaInput )
        {
          vevaciousPlusPlus.AppendResultsToLhaFile( runPointInput );
        }
      }
      else if( parameterElement->first == "ParameterPointSet" )
      {
        while( xmlParser.ReadNextElement() )
        {
          if( xmlParser.CurrentName() == "InputFolder" )
          {
            inputFolder = xmlParser.TrimmedCurrentBody();
          }
          else if( xmlParser.CurrentName() == "OutputFolder" )
          {
            outputFolder = xmlParser.TrimmedCurrentBody();
          }
          else if( xmlParser.CurrentName() == "AppendLhaOutputToLhaInput" )
          {
            appendLhaOutputToLhaInput = true;
          }
        }
        if( outputFolder.empty() )
        {
          std::cout
          << std::endl
          << "The <OutputFolder> content must not be an empty string! Please"
          << " use \"./\" for the current working folder.";
          std::cout << std::endl;

          return EXIT_FAILURE;
        }

        VevaciousPlusPlus::FilePlaceholderManager placeholderManager( "",
                                                                ".placeholder",
                                                                     ".vout" );
        placeholderManager.PrepareFilenames( inputFolder,
                                             outputFolder,
                                             outputFolder );

        while( placeholderManager.HoldNextPlace() )
        {
          vevaciousPlusPlus.RunPoint( placeholderManager.CurrentInput() );
          vevaciousPlusPlus.WriteResultsAsXmlFile(
                                          placeholderManager.CurrentOutput() );
          if( appendLhaOutputToLhaInput )
          {
            vevaciousPlusPlus.AppendResultsToLhaFile(
                                           placeholderManager.CurrentInput() );
          }
        }
      }
    }
  }

  std::cout
  << std::endl
  << "Vevacious finished running.";
  std::cout << std::endl;

  // this was a triumph! I'm making a note here:
  return EXIT_SUCCESS;
}
