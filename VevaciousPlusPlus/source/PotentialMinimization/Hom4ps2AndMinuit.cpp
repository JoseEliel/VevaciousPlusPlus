/*
 * Hom4ps2AndMinuit.cpp
 *
 *  Created on: Apr 7, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  Hom4ps2AndMinuit::Hom4ps2AndMinuit(
                      HomotopyContinuationReadyPolynomial& polynomialPotential,
                                      std::string const& pathToHom4ps2,
                                      std::string const homotopyType ) :
    HomotopyContinuationAndGradient( polynomialPotential ),
    polynomialPotential( polynomialPotential ),
    potentialForMinuit( polynomialPotential ),
    pathToHom4ps2( pathToHom4ps2 ),
    homotopyType( homotopyType ),
    variableNamer( 4,
                   '0',
                   6,
                   2,
                   "v" ),
    complexSolutions(),
    variableNames(),
    purelyRealSolutionSets(),
    minuitManager( potentialForMinuit )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "Hom4ps2AndMinuit::Hom4ps2AndMinuit( ..., \"" << pathToHom4ps2
    << "\", \"" << homotopyType << "\" )";
    std::cout << std::endl;/**/
  }

  Hom4ps2AndMinuit::~Hom4ps2AndMinuit()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "Hom4ps2AndMinuit::~Hom4ps2AndMinuit()";
    std::cout << std::endl;/**/
  }


  // This uses HOM4PS2 to fill purelyRealSolutionSets with all the extrema
  // of polynomialPotential.TargetPolynomialGradient().
  void Hom4ps2AndMinuit::FindTreeLevelExtrema()
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
    ParseHom4ps2Output( "./data.roots" );

    // Now we return to the original working directory so as to avoid confusing
    // the user.
    directoryChangeSuccess = chdir( originalWorkingDirectory );
    if( 0 != directoryChangeSuccess )
    {
      throw std::runtime_error(
                      "could not change directory back to initial directory" );
    }
  }

  // This should find the minimum at temperature temperatureInGev nearest to
  // minimumToAdjust (which is assumed to be a minimum of the potential at a
  // different temperature).
  PotentialMinimum Hom4ps2AndMinuit::AdjustMinimumForTemperature(
                                       PotentialMinimum const& minimumToAdjust,
                                                double const temperatureInGev )
  {
    potentialForMinuit.SetTemperature( temperatureInGev );
    return
    PotentialMinimum( minuitManager( minimumToAdjust.FieldConfiguration() ) );
  }

  void Hom4ps2AndMinuit::WriteHom4p2Input(
                                      std::string const& hom4ps2InputFilename )
  {
    polynomialPotential.PreparePolynomialHomotopyContinuation();
    std::vector< PolynomialSum > const&
    targetSystem( polynomialPotential.TargetPolynomialGradient() );

    variableNames.clear();
    for( unsigned int whichVariable( 0 );
         whichVariable < targetSystem.size();
         ++whichVariable )
    {
      variableNames.push_back( variableNamer.intToString( whichVariable ) );
      nameToIndexMap[ variableNames.back() ] = whichVariable;
    }

    std::ofstream hom4ps2Input( hom4ps2InputFilename.c_str() );
    hom4ps2Input << "{\n";
    for( std::vector< PolynomialSum >::const_iterator
         whichConstraint( targetSystem.begin() );
         whichConstraint < targetSystem.end();
         ++whichConstraint )
    {
      hom4ps2Input
      << whichConstraint->AsStringAtCurrentScale( variableNames ) << ";\n";
    }
    hom4ps2Input << "}\n";
    hom4ps2Input.close();
  }

  void Hom4ps2AndMinuit::ParseHom4ps2Output(
                                     std::string const& hom4ps2OutputFilename )
  {
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
        polynomialPotential.AppendPureRealSolutionAndValidSignFlips(
                                                         candidateRealSolution,
                                                      purelyRealSolutionSets,
                                                                     2.0 );
      }
    }

    // debugging:
    /**/std::cout << std::endl << "debugging:"
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
      std::cout << " } with tree-level depth "
      << polynomialPotential.QuickApproximation( *realSolution ) << std::endl;
    }
    std::cout << "}" << std::endl;
    std::cout << std::endl;/**/

    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "still need to run MINUIT!";
    std::cout << std::endl;/**/


    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "still need to move this out of ParseHom4ps2Output!";
    std::cout << std::endl;/**/

    MinuitManager minuitManager( potentialForMinuit );
    for( std::vector< std::vector< double > >::iterator
         realSolution( purelyRealSolutionSets.begin() );
         realSolution < purelyRealSolutionSets.end();
         ++realSolution )
    {
      /*ROOT::Minuit2::FunctionMinimum
      foundMinimum( minuitManager( *realSolution ) );*/

      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "foundMinimum = " << foundMinimum << std::endl;
      std::cout << std::endl;*/


      // placeholder:
      /**/std::cout << std::endl
      << "Placeholder: "
      << "still need to do something with ROOT::Minuit2::FunctionMinimum";
      std::cout << std::endl;/**/
    }
  }

  // This uses Minuit2 to minimize potentialForMinuit starting from the
  // values in purelyRealSolutionSets.
  void Hom4ps2AndMinuit::RollAndSortExtrema( )
  {
    double const
    thresholdSepartionSquared( ( 0.01 * dsbVacuum.LengthSquared() ) + 1.0 );
    double const thresholdSepartion( sqrt( thresholdSepartionSquared ) );
    PotentialMinimum foundMinimum;
    for( std::vector< std::vector< double > >::iterator
         realSolution( purelyRealSolutionSets.begin() );
         realSolution < purelyRealSolutionSets.end();
         ++realSolution )
    {
      foundMinimum = minuitManager( *realSolution );
      foundMinima.push_back( foundMinimum );

      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "found [" << foundMinimum.AsDebuggingString() << "]"
      << std::endl;
      std::cout << std::endl;/**/

      if( ( ( foundMinimum.FunctionValue() + foundMinimum.FunctionError() )
            < dsbVacuum.FunctionValue() )
          &&
          ( foundMinimum.SquareDistanceTo( dsbVacuum )
            > thresholdSepartionSquared )
          &&
          ( IsNotPhaseRotationOfDsbVacuum( foundMinimum,
                                           thresholdSepartion ) ) )
      {
        if( panicVacua.empty()
            ||
            ( foundMinimum.SquareDistanceTo( dsbVacuum )
              < panicVacuum.SquareDistanceTo( dsbVacuum ) ) )
        {
          panicVacuum = foundMinimum;
        }
        panicVacua.push_back( foundMinimum );
      }
    }

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "dsbVacuum = " << dsbVacuum.AsDebuggingString()
    << std::endl
    << "panicVacuum = " << panicVacuum.AsDebuggingString()
    << std::endl
    << "panicVacua.size() = " << panicVacua.size();
    for( unsigned int panicIndex( 0 );
         panicIndex < panicVacua.size();
         ++panicIndex )
    {
      std::cout << std::endl << "panicVacua[ " << panicIndex << " ] = "
      << panicVacua[ panicIndex ].AsDebuggingString();
    }
    std::cout << std::endl;
    std::cout << std::endl;/**/
  }

} /* namespace VevaciousPlusPlus */
