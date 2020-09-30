/*
 * VirtualSimpleLhaParser.hpp
 *  Created on: May 25 2018
 *      Author: José Eliel Camargo-Molina
 *  Derived class from simpleLhaparser,Extended to accept the creation of LhaBlockSet
 *  in code through the ReadBlock function  
 *
 */

#ifndef VIRTUALSIMPLELHAPARSER_HPP_
#define VIRTUALSIMPLELHAPARSER_HPP_

#include <vector>
#include <string>
#include <utility>
#include <cstddef>
#include <list>
#include <cctype>
#include <map>
#include "LHPC/Utilities/ParsingUtilities.hpp"
#include "LHPC/SimpleLhaParser.hpp"
#include <sstream>
#include <stdexcept>
#include <fstream>
#include "boost/format.hpp"


namespace VevaciousPlusPlus
{

  class VirtualSimpleLhaParser: public LHPC::SimpleLhaParser
    {
  public:
   
    // This opens the file with name fileName and parses it into blocks.
    void ReadFile( std::string const& fileName );
    
    
   // This prepares for reading in a parameter point as a set of blocks read with Read Block
    void DeleteBlocks();

	// This reads a block and creates the appropriate LhaBlock object storing that value. 
    // The Block name should be uppercase. The parameter values are given in a vector 
    // and a bool determines whether the slha indices should be sequential or matrix like 
 
    void ReadBlock(std::string const& uppercaseBlockName, double const scale, 
    			   std::vector<std::pair<int,double>> const& parameters, 
    			   int const dimension);
    
  };

  // This opens the file with name fileName and parses it into blocks.
  inline void VirtualSimpleLhaParser::ReadFile( std::string const& fileName )
  {
  	if(fileName != "internal" && fileName != "global" && fileName != "nearest")
  	{ 
     ResetForNewFile();
     std::string readLine( "" );
     std::ifstream fileStream( fileName.c_str() );
     if( !(fileStream.is_open()) )
     {
       std::stringstream errorBuilder;
       errorBuilder << "Could not open file named \"" << fileName << "\".";
       throw std::runtime_error( errorBuilder.str() );
     }
     while( std::getline( fileStream, readLine ) )
     {
       ParseLine( readLine );
     }
     fileStream.close();
    }
  }
  
// This prepares for reading in a parameter point as a set of blocks that have been read with Read Block

  inline void VirtualSimpleLhaParser::DeleteBlocks()
  {
   ResetForNewFile();
  }
  
// This reads a block and creates the appropriate LhaBlock object storing that value. 
// The Block name should be uppercase. The parameter values are given in a vector 
// and the dimension is given, with 1 is given the block is read as a list, if n is given
// the vector is read as a n x n matrix in sequential order row by row. 
 
 inline void VirtualSimpleLhaParser::ReadBlock(std::string const& uppercaseBlockName, 
 									 		   double const scale,
 									 		   std::vector<std::pair<int,double>> const& parameters, 
 									 		   int const dimension)
  { 
   currentBlockSet = BlockSetForName( uppercaseBlockName );
   currentBlock = AddNewBlockWithExplicitScale( scale );
   std::string linetopushback;
   if (dimension == 1)
   {
     for(unsigned int currentParameter=0;currentParameter < parameters.size(); currentParameter=currentParameter+1 ) 
    {
    linetopushback = (boost::format("%s %s") % parameters[currentParameter].first % parameters[currentParameter].second).str();  
    currentBlock->AddLine( linetopushback ); 
    }
   } 
   else
   {
   	 int row = 0;
     int column = 0 ;
     for(unsigned int currentParameter=0;currentParameter< parameters.size(); currentParameter=currentParameter+1 ) 
    {
 	row = parameters[currentParameter].first / 10; // Get first digit 
 	column= parameters[currentParameter].first % 10; // Get second digit
    linetopushback = (boost::format("%s %s %s") %  row % column % parameters[currentParameter].second).str();  
    currentBlock->AddLine( linetopushback ); 
    }
   }
  }
 
} /* namespace VevaciousPlusPlus */

#endif /* VIRTUALSIMPLELHAPARSER_HPP_ */
