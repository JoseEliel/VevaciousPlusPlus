/*
 * Hom4ps2Runner.cpp
 *
 *  Created on: Nov 23, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialMinimization/HomotopyContinuation/Hom4ps2Runner.hpp"

namespace VevaciousPlusPlus
{
  std::string const Hom4ps2Runner::fieldNamePrefix( "fv" );

  Hom4ps2Runner::Hom4ps2Runner( std::string const& pathToHom4ps2,
                                std::string const& homotopyType,
                                double const resolutionSize ) :
    pathToHom4ps2( pathToHom4ps2 ),
    homotopyType( homotopyType ),
    resolutionSize( resolutionSize )
  {
    // This constructor is just an initialization list.
  }

  Hom4ps2Runner::~Hom4ps2Runner()
  {
    // This does nothing.
  }


  // This uses HOM4PS2 to fill startingPoints with all the extrema of
  // targetSystem.TargetPolynomialGradient().
  void Hom4ps2Runner::operator()(
                      std::vector< PolynomialConstraint > const& systemToSolve,
                  std::vector< std::vector< double > >& systemSolutions ) const
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
    std::vector< std::string > variableNames( systemToSolve.size(),
                                              "" );
    std::map< std::string, size_t > nameToIndexMap;
    WriteHom4p2Input( systemToSolve,
                      variableNames,
                      nameToIndexMap,
                      hom4ps2InputFilename );

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
                        systemSolutions,
                        variableNames,
                        nameToIndexMap,
                        systemToSolve );

    // Now we return to the original working directory so as to avoid confusing
    // the user.
    directoryChangeSuccess = chdir( originalWorkingDirectory );
    if( 0 != directoryChangeSuccess )
    {
      throw std::runtime_error(
                      "could not change directory back to initial directory" );
    }
  }

  // This sets up the variable names in variableNames and nameToIndexMap,
  // then writes systemToSolve using these names in the correct form for
  // HOM4PS2 in a file with name hom4ps2InputFilename.
  void Hom4ps2Runner::WriteHom4p2Input(
                      std::vector< PolynomialConstraint > const& systemToSolve,
                                     std::vector< std::string >& variableNames,
                               std::map< std::string, size_t >& nameToIndexMap,
                                std::string const& hom4ps2InputFilename ) const
  {
    size_t const numberOfFields( systemToSolve.size() );
    variableNames.resize( numberOfFields );
    std::stringstream nameBuilder;
    nameBuilder << numberOfFields;
    size_t const numberOfDigits( nameBuilder.str().size() );
    nameBuilder.fill( '0' );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      nameBuilder.str( "" );
      nameBuilder.width( numberOfDigits );
      nameBuilder << ( fieldIndex + 1 );
      variableNames[ fieldIndex ] = ( fieldNamePrefix + nameBuilder.str() );
      nameToIndexMap[ variableNames[ fieldIndex ] ] = fieldIndex;
    }

    std::ofstream hom4ps2Input( hom4ps2InputFilename.c_str() );
    hom4ps2Input << "{\n";
    for( std::vector< PolynomialConstraint >::const_iterator
         constraintToWrite( systemToSolve.begin() );
         constraintToWrite != systemToSolve.end();
         ++constraintToWrite )
    {
      hom4ps2Input
      << WriteConstraint( *constraintToWrite,
                          variableNames ) << ";\n";
    }
    hom4ps2Input << "}\n";
    hom4ps2Input.close();
  }

  // This returns the constraint as a string of terms joined by '+' or '-'
  // appropriately, where each term is of the form
  // coefficient " * " variableName[ fieldIndex ] "^" appropriate power
  // (without writing any power part if the power is only 1, and without
  // writing the field name at all if its power is 0).
  std::string Hom4ps2Runner::WriteConstraint(
                                 PolynomialConstraint const& constraintToWrite,
                        std::vector< std::string > const& variableNames ) const
  {
    std::stringstream stringBuilder;
    bool firstTermWritten( false );
    for( std::vector< FactorWithPowers >::const_iterator
         factorWithPowers( constraintToWrite.begin() );
         factorWithPowers != constraintToWrite.end();
         ++factorWithPowers )
    {
      if( factorWithPowers->first != 0.0 )
      {
        if( !firstTermWritten )
        {
          stringBuilder << factorWithPowers->first;
        }
        else if( factorWithPowers->first < 0.0 )
        {
          stringBuilder << " - " << -(factorWithPowers->first);
        }
        else
        {
          stringBuilder << " + " << factorWithPowers->first;
        }

        for( size_t fieldIndex( 0 );
             fieldIndex < factorWithPowers->second.size();
             ++fieldIndex )
        {
          if( factorWithPowers->second[ fieldIndex ] > 0 )
          {
            stringBuilder << " * " << variableNames[ fieldIndex ];
            if( factorWithPowers->second[ fieldIndex ] > 1 )
            {
              stringBuilder << "^" << factorWithPowers->second[ fieldIndex ];
            }
          }
        }
        firstTermWritten = true;
      }
    }
    return stringBuilder.str();
  }

  void
  Hom4ps2Runner::ParseHom4ps2Output( std::string const& hom4ps2OutputFilename,
                  std::vector< std::vector< double > >& purelyRealSolutionSets,
                               std::vector< std::string > const& variableNames,
                         std::map< std::string, size_t > const& nameToIndexMap,
               std::vector< PolynomialConstraint > const& systemToSolve ) const
  {
    std::cout
    << std::endl
    << "-----------------" << std::endl << "Parsing output from HOM4PS2."
    << std::endl ;
    std::cout << std::endl;

    std::vector< std::complex< long double > > complexSolutions;
    BOL::CommentedTextParser tadpoleSolutionsFile( "###",
                                                   false );
    bool successfulOperation( tadpoleSolutionsFile.openFile(
                                                     hom4ps2OutputFilename ) );
    if( !successfulOperation )
    {
      std::stringstream errorBuilder;
      errorBuilder
      << "Could not open " <<  hom4ps2OutputFilename
      << " to write input file for HOM4PS2.";
      throw std::runtime_error( errorBuilder.str() );
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
      if( lineString[ 0 ] == '(' )
      {
        // If the above condition is satisfied, lineString now contains the
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
      else if( lineString == "The order of variables :" )
      {
        break;
      }
    }

    size_t const numberOfVariables( variableNames.size() );
    size_t const
    numberOfComplexSolutions( complexSolutions.size() / numberOfVariables );
    if( ( numberOfComplexSolutions * numberOfVariables )
        != complexSolutions.size() )
    {
      std::stringstream errorBuilder;
      errorBuilder << "Parsed " << complexSolutions.size()
      << " values from the output of HOM4PS2, but it should have been a"
      << " multiple of " << numberOfVariables << " as there are "
      << numberOfVariables << " variables.";
      throw std::runtime_error( errorBuilder.str() );
    }

    // At this point, the line "The order of variables :" should have been
    // found. If it hasn't, the file is malformed, but we carry on regardless,
    // looking for the variables in order:
    std::vector< size_t > indexOrder( numberOfVariables,
                                      -1 );
    std::map< std::string, size_t >::const_iterator indexFinder;
    size_t variableIndex( 0 );
    while( tadpoleSolutionsFile.readNextNonEmptyLineOfFileWithoutComment(
                                                                 lineString ) )
    {
      if( lineString == "===============>   HOM4PS-2.0   <===============" )
      {
        break;
      }
      indexFinder = nameToIndexMap.find(
                        LHPC::ParsingUtilities::TrimWhitespaceFromFrontAndBack(
                                                                lineString ) );
      if( indexFinder == nameToIndexMap.end() )
      {
        std::stringstream errorBuilder;
        errorBuilder
        << "Hom4ps2Runner::ParseHom4ps2Output managed to find a variable name"
        << " (\""
        << LHPC::ParsingUtilities::TrimWhitespaceFromFrontAndBack( lineString )
        << "\") in the HOM4PS2 output which it did not create.";
        throw
        std::runtime_error( errorBuilder.str() );
      }
      indexOrder[ variableIndex ] = indexFinder->second;
      ++variableIndex;
    }

    // Now we divide complexSolutions into sets of numberOfVariables complex
    // values, and any which are purely real (within a tolerance of
    // resolutionSize) are kept, along with any valid sign-flip variations.
    std::vector< double > candidateRealSolution( numberOfVariables,
                                                 0.0 );
    bool solutionIsReal( true );
    size_t currentStartIndex( 0 );
    for( size_t solutionIndex( 0 );
         solutionIndex < numberOfComplexSolutions;
         ++solutionIndex )
    {
      solutionIsReal = true;
      for( size_t variableIndex( 0 );
           variableIndex < numberOfVariables;
           ++variableIndex )
      {
        if( abs( complexSolutions[ currentStartIndex + variableIndex ].imag() )
            > resolutionSize )
        {
          solutionIsReal = false;
          break;
        }
        else
        {
          candidateRealSolution[ indexOrder[ variableIndex ] ]
          = complexSolutions[ currentStartIndex + variableIndex ].real();
        }
      }
      currentStartIndex += numberOfVariables;
      if( solutionIsReal )
      {
        AppendSolutionAndValidSignFlips( candidateRealSolution,
                                         purelyRealSolutionSets,
                                         systemToSolve,
                                         resolutionSize );
      }
    }

    std::cout
    << std::endl
    << "-----------------" << std::endl << "Parsed "
    << ( complexSolutions.size() / numberOfVariables )
    << " complex solutions from HOM4PS2. After trying sign-flip variations,"
    << " returning " << purelyRealSolutionSets.size()
    << " purely real solutions."
    << std::endl;
    std::cout << std::endl;
  }

} /* namespace VevaciousPlusPlus */
