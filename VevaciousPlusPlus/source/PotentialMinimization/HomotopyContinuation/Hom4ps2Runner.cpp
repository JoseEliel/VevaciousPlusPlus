/*
 * Hom4ps2Runner.cpp
 *
 *  Created on: Apr 14, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  Hom4ps2Runner::Hom4ps2Runner( PolynomialGradientTargetSystem& targetSystem,
                                std::string const& pathToHom4ps2,
                                std::string const homotopyType ) :
    HomotopyContinuationSolver( targetSystem ),
    targetSystem( targetSystem ),
    pathToHom4ps2( pathToHom4ps2 ),
    homotopyType( homotopyType ),
    variableNamer( 4,
                   '0',
                   6,
                   2,
                   "v" ),
    complexSolutions(),
    variableNames(),
    nameToIndexMap()
  {
    // This constructor is just an initialization list.
  }

  Hom4ps2Runner::~Hom4ps2Runner()
  {
    // This does nothing.
  }


  // This uses HOM4PS2 to fill purelyRealSolutionSets with all the extrema of
  // targetSystem.TargetPolynomialGradient().
  void Hom4ps2Runner::FindTreeLevelExtrema(
                 std::vector< std::vector< double > >& purelyRealSolutionSets )
  {
    char originalWorkingDirectory[ PATH_MAX ];
    if( NULL == getcwd( originalWorkingDirectory,
                        PATH_MAX ) )
    {
      std::cout
      << std::endl
      << "Error! unable to determine current working directory! (necessary,"
      << " since this program needs to change directory to the directory where"
      << " the hom4ps2 executable is, since unfortunately HOM4PS2 runs with"
      << " relative paths; this program returns to where it was called though,"
      << " to make batch calls easier.)";
      std::cout << std::endl;
      throw std::runtime_error(
                             "could not determine current working directory" );
    }
    int directoryChangeSuccess( chdir( pathToHom4ps2.c_str() ) );
    if( 0 != directoryChangeSuccess )
    {
      throw std::runtime_error(
                           "could not change directory to HOM4PS2 directory" );
    }
    std::string hom4ps2InputFilename( "./VevaciousHomotopyContinuation.txt" );

    WriteHom4p2Input( hom4ps2InputFilename );

    std::cout
    << std::endl
    << "Running HOM4PS2!" << std::endl << "-----------------" << std::endl;
    std::cout << std::endl;

    std::string systemCommand( "rm ./bin/input.num" );
    BOL::UsefulStuff::runSystemCommand( systemCommand );
    systemCommand.assign( "/bin/bash -c \"./hom4ps2 " );
    systemCommand.append( hom4ps2InputFilename );
    systemCommand.append(  " <<< " );
    systemCommand.append( homotopyType );
    systemCommand.append( "\"" );
    BOL::UsefulStuff::runSystemCommand( systemCommand );

    // At this point, we are in the directory with hom4ps2 & data.roots, so
    // now we fill purelyRealSolutionSets.
    ParseHom4ps2Output( "./data.roots",
                        purelyRealSolutionSets );

    // Now we return to the original working directory so as to avoid confusing
    // the user.
    directoryChangeSuccess = chdir( originalWorkingDirectory );
    if( 0 != directoryChangeSuccess )
    {
      throw std::runtime_error(
                      "could not change directory back to initial directory" );
    }
  }


  void
  Hom4ps2Runner::WriteHom4p2Input( std::string const& hom4ps2InputFilename )
  {
    variableNames.clear();
    for( unsigned int whichVariable( 0 );
         whichVariable < targetSystem.TargetSystem().size();
         ++whichVariable )
    {
      variableNames.push_back( variableNamer.intToString( whichVariable ) );
      nameToIndexMap[ variableNames.back() ] = whichVariable;
    }

    std::ofstream hom4ps2Input( hom4ps2InputFilename.c_str() );
    hom4ps2Input << "{\n";
    for( std::vector< PolynomialSum >::const_iterator
         whichConstraint( targetSystem.TargetSystem().begin() );
         whichConstraint < targetSystem.TargetSystem().end();
         ++whichConstraint )
    {
      hom4ps2Input
      << whichConstraint->AsStringAtCurrentScale( variableNames ) << ";\n";
    }
    hom4ps2Input << "}\n";
    hom4ps2Input.close();
  }

  void
  Hom4ps2Runner::ParseHom4ps2Output( std::string const& hom4ps2OutputFilename,
                 std::vector< std::vector< double > >& purelyRealSolutionSets )
  {
    std::cout
    << std::endl
    << "-----------------" << std::endl << "Parsing output from HOM4PS2."
    << std::endl ;
    std::cout << std::endl;

    complexSolutions.clear();
    purelyRealSolutionSets.clear();
    BOL::CommentedTextParser tadpoleSolutionsFile( "###",
                                                   false );
    bool successfulOperation( tadpoleSolutionsFile.openFile(
                                                     hom4ps2OutputFilename ) );
    if( !successfulOperation )
    {
      throw std::runtime_error( "could not open " + hom4ps2OutputFilename );
    }
    // First we pick out lines corresponding to solutions until we find
    // "The order of variables :" which comes after all solutions have been
    // printed.
    std::string lineString( "" );
    std::complex< long double > currentComplexNumber( 0.0L,
                                                      0.0L );
    std::stringstream parsingStream;
    while( tadpoleSolutionsFile.readNextNonEmptyLineOfFileWithoutComment(
                                                                 lineString ) )
    {
      if( '(' == lineString[ 0 ] )
      {
        // If the above conditions are satisfied, lineString now contains the
        // root as a complex number, in the form where zero is
        // "(  0.0000000000000000E+00 ,  0.0000000000000000E+00)"
        BOL::StringParser::substituteCharacterWith( lineString,
                                                    ',',
                                                    ' ' );
        parsingStream.clear();
        parsingStream.str( lineString.substr( 1,
                                              ( lineString.size() - 2 ) ) );
        parsingStream
        >> currentComplexNumber.real() >> currentComplexNumber.imag();
        complexSolutions.push_back( currentComplexNumber );
      }
      else if( 0 == lineString.compare( "The order of variables :" ) )
      {
        break;
      }
    }
    // At this point, the line "The order of variables :" should have been
    // found. If it hasn't, the file is malformed, but we carry on regardless,
    // looking for the variables in order:
    unsigned int numberOfVariables( variableNames.size() );
    std::vector< unsigned int > indexOrder( numberOfVariables );
    unsigned int whichVariable( 0 );
    while( tadpoleSolutionsFile.readNextNonEmptyLineOfFileWithoutComment(
                                                                 lineString ) )
    {
      if( 0 == lineString.compare(
                         "===============>   HOM4PS-2.0   <===============" ) )
      {
        break;
      }
      indexOrder[ whichVariable ]
      = nameToIndexMap[ BOL::StringParser::trimFromFrontAndBack(
                                                                lineString ) ];
      ++whichVariable;
    }

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "indexOrder = ";
    for( std::vector< unsigned int >::iterator
         whichIndex( indexOrder.begin() );
         whichIndex < indexOrder.end();
         ++whichIndex )
    {
      std::cout << std::endl << *whichIndex;
    }
    std::cout << std::endl
    << "complexSolutions = { " << std::endl;
    for( std::vector< std::complex< long double > >::iterator
         complexValue( complexSolutions.begin() );
         complexValue < complexSolutions.end();
         ++complexValue )
    {
      std::cout << " ( " << complexValue->real() << ", "
      << complexValue->imag() << " i )";
    }
    std::cout << "}" << std::endl;
    std::cout << std::endl;*/

    std::vector< std::complex< double > >
    candidateRealSolution( numberOfVariables );
    unsigned int solutionIndex;
    for( unsigned int complexIndex( 0 );
         complexIndex < complexSolutions.size();
         ++complexIndex )
    {
      solutionIndex = ( complexIndex % numberOfVariables );
      candidateRealSolution[ indexOrder[ solutionIndex ] ].real()
      = (double)(complexSolutions[ complexIndex ].real());
      candidateRealSolution[ indexOrder[ solutionIndex ] ].imag()
      = (double)(complexSolutions[ complexIndex ].imag());
      if( solutionIndex == ( numberOfVariables - 1 ) )
      {
        targetSystem.AppendPureRealSolutionAndValidSignFlips(
                                                         candidateRealSolution,
                                                      purelyRealSolutionSets,
                                                              1.0 );
      }
    }

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "purelyRealSolutionSets = {" << std::endl;
    for( std::vector< std::vector< double > >::iterator
         realSolution( purelyRealSolutionSets.begin() );
         realSolution < purelyRealSolutionSets.end();
         ++realSolution )
    {
      std::cout << "{";
      for( std::vector< double >::iterator
           realValue( realSolution->begin() );
           realValue < realSolution->end();
           ++realValue )
      {
        std::cout << " " << *realValue;
      }
      std::cout << " }" << std::endl;
    }
    std::cout << "}" << std::endl;
    std::cout << std::endl;*/
  }

} /* namespace VevaciousPlusPlus */
