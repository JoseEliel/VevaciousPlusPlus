/*
 * Hom4ps2Runner.cpp
 *
 *  Created on: Nov 23, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialMinimization/HomotopyContinuation/Hom4ps2Runner.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/random_generator.hpp>

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
    if( !(( homotopyType == "1" )
          ||
          ( homotopyType == "2" )) )
    {
      std::stringstream errorBuilder;
      errorBuilder << "HOM4PS2 argument must be either"
      << " \"1\" (for \"The polyhedral homotopy\")"
      << " or \"2\" (for \"The classical linear homotopy\")";
      throw std::runtime_error( errorBuilder.str() );
    }
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
    // Here we change to HOM4PS2 directory
    int directoryChangeSuccess( chdir( pathToHom4ps2.c_str() ) );
    if( 0 != directoryChangeSuccess )
    {
      throw std::runtime_error(
                           "could not change directory to HOM4PS2 directory" );
    }

    // Here we create a unique name for a folder within HOM4PS2's folder to run
    // the point.

    std::string pathname = boost::lexical_cast<std::string>(boost::uuids::random_generator()());
    std::string absolutepathname = pathToHom4ps2 + "/" + pathname;

    // Here we make to the unique directory

    std::string systemCommand( "mkdir "+ pathname );

    int systemReturn( system( systemCommand.c_str() ) );
    if( systemReturn == -1 )
    {
      std::stringstream errorBuilder;
      errorBuilder << "System could not execute \"" << systemCommand << "\".";
      throw std::runtime_error( errorBuilder.str() );
    }

    // Here we go inside the unique directory

    int tempdirectoryChangeSuccess( chdir( absolutepathname.c_str() ) );
    if( 0 != tempdirectoryChangeSuccess )
    {
      throw std::runtime_error(
              "could not change directory to the temporary folder within HOM4PS2 directory" );
    }

    // Here we make a symlink to the folder with HOM4PS2 binary files (bin/) within the unique
    // directory, otherwise HOM4PS2 main executable won't work as it always searchers for the
    // binary files in ./bin (relative path)

    systemCommand.assign( "ln -s ../bin ./bin" );

    systemReturn = system( systemCommand.c_str() ) ;
    if( systemReturn == -1 )
    {
      std::stringstream errorBuilder;
      errorBuilder << "System could not execute \"" << systemCommand << "\".";
      throw std::runtime_error( errorBuilder.str() );
    }

    std::string hom4ps2InputFilename("VevaciousHomotopyContinuation.txt" );


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

    // This is to avoid a bug in HOM4PS2 where it tries to run this file instead of the input.

    systemCommand.assign( "rm ../bin/input.num" );
    systemReturn = system( systemCommand.c_str() );
    if( systemReturn == -1 )
    {
      std::stringstream errorBuilder;
      errorBuilder << "System could not execute \"" << systemCommand << "\".";
      throw std::runtime_error( errorBuilder.str() );
    }
    // Running HOM4PS2
    systemCommand.assign( "/bin/bash -c \"../hom4ps2 " );
    systemCommand.append( hom4ps2InputFilename );
    systemCommand.append(  " <<< " );
    systemCommand.append( homotopyType );
    systemCommand.append( "\"" );
    systemReturn = system( systemCommand.c_str() );
    if( systemReturn == -1 )
    {
      std::stringstream errorBuilder;
      errorBuilder << "System could not execute \"" << systemCommand << "\".";
      throw std::runtime_error( errorBuilder.str() );
    }

    // At this point, we are in the unique directory with hom4ps2 & data.roots
    // one directory above, so now we fill purelyRealSolutionSets.

    ParseHom4ps2Output( "./data.roots",
                        systemSolutions,
                        variableNames,
                        nameToIndexMap,
                        systemToSolve );

    systemCommand.assign( "rm ../bin/input.num" );
    systemReturn = system( systemCommand.c_str() );
    if( systemReturn == -1 )
    {
      std::stringstream errorBuilder;
      errorBuilder << "System could not execute \"" << systemCommand << "\".";
      throw std::runtime_error( errorBuilder.str() );
    }
    // Now we return to the original working directory so as to avoid confusing
    // the user.
    directoryChangeSuccess = chdir( originalWorkingDirectory );
    if( 0 != directoryChangeSuccess )
    {
      throw std::runtime_error(
                      "could not change directory back to initial directory" );
    }

    systemCommand.assign( "rm -rf " + absolutepathname );

    systemReturn = system( systemCommand.c_str() ) ;
    if( systemReturn == -1 )
    {
      std::stringstream errorBuilder;
      errorBuilder << "System could not execute \"" << systemCommand << "\".";
      throw std::runtime_error( errorBuilder.str() );
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
      hom4ps2Input << WriteConstraint( *constraintToWrite,
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
    std::ifstream tadpoleSolutionsFile( hom4ps2OutputFilename.c_str() );
    if( !(tadpoleSolutionsFile.good()) )
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
    while( std::getline( tadpoleSolutionsFile,
                         lineString ).good() )
    {
      if( lineString[ 0 ] == '(' )
      {
        // If the above condition is satisfied, lineString now contains the
        // root as a complex number, in the form where zero is
        // "(  0.0000000000000000E+00 ,  0.0000000000000000E+00)"
		double Re,Im;
        LHPC::ParsingUtilities::ReplaceAllCharacters( lineString,
                                                      ',',
                                                      ' ' );
        parsingStream.clear();
        parsingStream.str( lineString.substr( 1,
                                              ( lineString.size() - 2 ) ) );
        parsingStream
        >> Re >> Im;
		currentComplexNumber.real(Re);
		currentComplexNumber.imag(Im);
        complexSolutions.push_back( currentComplexNumber );
      }
      else if( LHPC::ParsingUtilities::TrimWhitespaceFromFrontAndBack(
                                                                   lineString )
               == "The order of variables :" )
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
      << ( ( complexSolutions.size() == 1 ) ? " value" : " values" )
      << " from the output of HOM4PS2, but it should have been a multiple of "
      << numberOfVariables << " as there are " << numberOfVariables
      << " variables.";
      throw std::runtime_error( errorBuilder.str() );
    }

    // At this point, the line "The order of variables :" should have been
    // found. If it hasn't, the file is malformed, but we carry on regardless,
    // looking for the variables in order:
    std::vector< size_t > indexOrder( numberOfVariables,
                                      -1 );
    std::map< std::string, size_t >::const_iterator indexFinder;
    size_t variableIndex( 0 );
    while( std::getline( tadpoleSolutionsFile,
                         lineString ).good() )
    {
      lineString
      = LHPC::ParsingUtilities::TrimWhitespaceFromFrontAndBack( lineString );
      if( lineString.empty() )
      {
        continue;
      }
      else if( lineString
               == "===============>   HOM4PS-2.0   <===============" )
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
        if( fabs(
                 complexSolutions[ currentStartIndex + variableIndex ].imag() )
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

    unsigned int const numberOfParsedComplexSolutions( complexSolutions.size()
                                                       / numberOfVariables );
    std::cout
    << std::endl
    << "-----------------" << std::endl << "Parsed "
    << numberOfParsedComplexSolutions
    << " complex solution"
    << ( ( numberOfParsedComplexSolutions == 1 ) ? "" : "s" )
    << " from HOM4PS2. After trying sign-flip variations,"
    << " returning " << purelyRealSolutionSets.size()
    << " purely real solution"
    << ( ( purelyRealSolutionSets.size() == 1 ) ? "." : "s." )
    << std::endl;
    std::cout << std::endl;
  }

} /* namespace VevaciousPlusPlus */
