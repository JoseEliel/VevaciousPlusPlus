/*
 * PHCRunner.cpp
 *
 *  Created on: Nov 22, 2017
 *      Author: Simon Geisler (simon.geisler94@gmail.com)
 */

#include "PotentialMinimization/HomotopyContinuation/PHCRunner.hpp"

namespace VevaciousPlusPlus
{
  std::string const PHCRunner::fieldNamePrefix( "fv" );

  PHCRunner::PHCRunner( std::string const& pathToPHC,
                                double const resolutionSize, unsigned const int taskcount ) :
    pathToPHC( pathToPHC ),
    resolutionSize( resolutionSize ),
    taskcount(  taskcount  )
  {
  }

  PHCRunner::~PHCRunner()
  {
    // This does nothing.
  }


  // This uses PHC to fill startingPoints with all the extrema of
  // targetSystem.TargetPolynomialGradient().
  void PHCRunner::operator()(
                      std::vector< PolynomialConstraint > const& systemToSolve,
                  std::vector< std::vector< double > >& systemSolutions ) const
  {

    std::string PHCInputFileName = pathToPHC + "VevaciousHomotopyContinuation.txt";
    std::string PHCOutputFilename = pathToPHC + "PHCOutput";
    std::vector< std::string > variableNames( systemToSolve.size(),
                                              "" );
    std::map< std::string, size_t > nameToIndexMap;
    WritePHCInput( systemToSolve,
                      variableNames,
                      nameToIndexMap,
                      PHCInputFileName );

    std::cout
    << std::endl
    << "Running PHC!" << std::endl << "-----------------" << std::endl;
        std::cout << std::endl;
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	int systemReturn(0);
	std::string systemCommand = "rm " + PHCOutputFilename;
	
	struct stat buffer; //Checking if file exists, fastest method.
	if(stat(PHCOutputFilename.c_str(), &buffer)==0){		
		systemReturn = system(systemCommand.c_str());
		if( systemReturn == -1 )
		{
		  std::stringstream errorBuilder;
		  errorBuilder << "System could not remove PHCOutputfile with \"" << systemCommand << "\".";
		  throw std::runtime_error( errorBuilder.str() );
		}
	}
	systemCommand.assign(pathToPHC + "phc -b -t" + std::to_string(taskcount)+ " "); //calls the blackbox solver
    systemCommand.append( PHCInputFileName );
	systemCommand.append(" ");
	systemCommand.append( PHCOutputFilename );
    systemReturn = system( systemCommand.c_str() );
    if( systemReturn == -1 )
    {
      std::stringstream errorBuilder;
      errorBuilder << "System could not run PHC with \"" << systemCommand << "\".";
      throw std::runtime_error( errorBuilder.str() );
    }
	std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now(); // we want to measure the elapsed time
	std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " ms" <<std::endl;
    // now we fill purelyRealSolutionSets.
	begin = std::chrono::steady_clock::now();
    ParsePHCOutput( PHCInputFileName,
                        systemSolutions,
                        variableNames,
                        nameToIndexMap,
                        systemToSolve );
	end= std::chrono::steady_clock::now();
	std::cout << "Parsing time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " ms" <<std::endl << std::endl  << "-----------------" << std::endl;
	
  }

  // This sets up the variable names in variableNames and nameToIndexMap,
  // then writes systemToSolve using these names in the correct form for
  // PHC in a file with name PHCInputFilename.
  void PHCRunner::WritePHCInput(
                      std::vector< PolynomialConstraint > const& systemToSolve,
                                     std::vector< std::string >& variableNames,
                               std::map< std::string, size_t >& nameToIndexMap,
                                std::string const& PHCInputFilename ) const
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

    std::ofstream PHCInput( PHCInputFilename.c_str() );
	
    PHCInput << numberOfFields << " " << systemToSolve.size() << std::endl;
    for( std::vector< PolynomialConstraint >::const_iterator
         constraintToWrite( systemToSolve.begin() );
         constraintToWrite != systemToSolve.end();
         ++constraintToWrite )
    {
      PHCInput << WritePHCConstraint( *constraintToWrite,
                                       variableNames ) << "\n";
    }
    PHCInput << "\n";
    PHCInput.close();
  }

  // This returns the constraint as a string of terms joined by '+' or '-'
  // appropriately, where each term is of the form
  // coefficient " * " variableName[ fieldIndex ] "^" appropriate power
  // (without writing any power part if the power is only 1, and without
  // writing the field name at all if its power is 0).
    std::string PHCRunner::WritePHCConstraint(
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
	stringBuilder << ";";
    return stringBuilder.str();;
  }

  void
  PHCRunner::ParsePHCOutput( std::string const& PHCInputFileName,
                  std::vector< std::vector< double > >& purelyRealSolutionSets,
                               std::vector< std::string > const& variableNames,
                         std::map< std::string, size_t > const& nameToIndexMap,
               std::vector< PolynomialConstraint > const& systemToSolve ) const
  {

    std::cout
    << std::endl
    << "-----------------" << std::endl << std::endl << "Parsing the solutions of PHCpack"
    << std::endl;
	
	size_t const numberOfVariables( variableNames.size() );
	//Reading from file
	std::ifstream t(PHCInputFileName); //PHC appends the final solutions to the Inputfile. The Outputfile contains further information, which isn't needed.
	std::string container((std::istreambuf_iterator<char>(t)),
							std::istreambuf_iterator<char>());
	//-----------------
	
	//Parsing container
	std::map<int,std::vector<double>, std::less<int>> solmap;
	//This is an optimized algorithm for big solution containers. It will save memory and running time.
	for(auto it = nameToIndexMap.begin(); it!= nameToIndexMap.end(); it++) //running over all fields
		{
		std::string doublepattern (it->first); //Fieldvalue names.
		doublepattern += "\\s+:\\s+"; //whitespaces, colon, whitespaces
		doublepattern += "(-?[0-9]+.[0-9]+E[+-][0-9]+)\\s+(-?[0-9]+.[0-9]+E[+-][0-9]+)"; //match[1]: Re in scientific double; whitespaces; match[2] : Im in scientific double
		std::regex pattern(doublepattern);
		std::sregex_iterator next(container.begin(), container.end(), pattern);
		std::sregex_iterator end;
		double Re,Im;
		int step(0);
		if(it == nameToIndexMap.begin()){ //With the first iteration we need to fill solmap with all Real occurences of the starting Variable
			while (next != end) { //Going through all matches.
				std::smatch match = *next;
				Re = std::stod(match[1]);
				Im = std::stod(match[2]);
				if(fabs(Im) < resolutionSize) solmap[step].push_back(Re);
				++step;
				++next;
					}
			}
		else {//Now we just check the matches which had real occurences before
			for(auto itm=solmap.begin(); itm != solmap.end();)
			{
				for(int k=0;k<(itm->first - step);k++) ++next;  //We don't want to iterate all over from the beginning every time, we just go through all entrys of solmap in one cumulative iteration
				step = itm->first;
				std::smatch match = *next;
				Re = std::stod(match[1]);
				Im = std::stod(match[2]);
				if(fabs(Im) < resolutionSize) {
					solmap[itm->first].push_back(Re);
					itm++;
				}
				else {
					auto itm2 = itm;
					itm++;
					solmap.erase(itm2);
				} //We want all fieldvalues to be real, if one isn't, then the vector is erased. We need a temp iterator to not mess up the map order.
			}
		}
	}
	//-------------------
	//Appending Solutions

	if(!(solmap.empty())){
		for(auto it = solmap.begin(); it !=solmap.end(); it++)
		{
			if((it->second).size() == numberOfVariables) {
				AppendSolutionAndValidSignFlips(it->second,
											 purelyRealSolutionSets,
											 systemToSolve,
											resolutionSize); //Sign flips, because why not.
			}
			else 
			{
				std::stringstream errorBuilder;
				errorBuilder << "There seems to be an error, while parsing the real solutions. Check the Output of PHCpack, maybe it's empty or faulty because of an error." << std::endl;
				throw std::runtime_error( errorBuilder.str() );
			}
		}
	}
	else {
		std::stringstream errorBuilder;
		errorBuilder << "No real solutions have been found. Check on your ResolutionSize or your input system." << std::endl;
		throw std::runtime_error( errorBuilder.str() );
		}
	
    unsigned int const numberOfParsedRealSolutions(solmap.size());
    std::cout
    << "Parsed "
    << numberOfParsedRealSolutions
    << " real solution"
    << ( ( numberOfParsedRealSolutions == 1 ) ? "" : "s" )
    << " from PHC. "<<std::endl << "After trying sign-flip variations,"
    << " returning " << purelyRealSolutionSets.size()
    << " purely real solution"
    << ( ( purelyRealSolutionSets.size() == 1 ) ? "." : "s." )
    << std::endl  << std::endl;
  }

} /* namespace VevaciousPlusPlus */
