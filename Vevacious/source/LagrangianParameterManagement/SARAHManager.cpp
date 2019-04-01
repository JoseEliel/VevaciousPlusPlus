/*
 * SARAHManager.cpp
 *
 *  Created on: Nov 3, 2015
 *      Author: Simon Geisler 
 */


#include "LagrangianParameterManagement/SARAHManager.hpp"

namespace VevaciousPlusPlus
{

  SARAHManager::SARAHManager(
                                          std::string const& validBlocksString,
                                           std::string const& minimumScaleType,
                                       std::string const& minimumScaleArgument,
                                             std::string const& fixedScaleType,
                                         std::string const& fixedScaleArgument,
                                           std::string const& maximumScaleType,
                                    std::string const& maximumScaleArgument ) :
    LesHouchesAccordBlockEntryManager( validBlocksString,
                                       minimumScaleType,
                                       minimumScaleArgument,
                                       fixedScaleType,
                                       fixedScaleArgument,
                                       maximumScaleType,
                                       maximumScaleArgument ),
    activeDerivedParameters()
  {
    RegisterDerivedParameters(derivedparameters);
  }

  SARAHManager::SARAHManager(
                                 std::set< std::string > const& validBlocksSet,
                                           std::string const& minimumScaleType,
                                       std::string const& minimumScaleArgument,
                                             std::string const& fixedScaleType,
                                         std::string const& fixedScaleArgument,
                                           std::string const& maximumScaleType,
                                    std::string const& maximumScaleArgument ) :
    LesHouchesAccordBlockEntryManager( validBlocksSet,
                                       minimumScaleType,
                                       minimumScaleArgument,
                                       fixedScaleType,
                                       fixedScaleArgument,
                                       maximumScaleType,
                                       maximumScaleArgument ),
    activeDerivedParameters()
  {
	  RegisterDerivedParameters(derivedparameters);
  }

  SARAHManager::SARAHManager( std::string const& xmlFileName ) :
    LesHouchesAccordBlockEntryManager( xmlFileName ),
    activeDerivedParameters()
  {
	  RegisterDerivedParameters(derivedparameters);
  }

  SARAHManager::~SARAHManager()
  {

    for( size_t deletionIndex( 0 );
         deletionIndex < activeDerivedParameters.size();
         ++deletionIndex )
    {
      delete activeDerivedParameters[ deletionIndex ];
    }
  }


  //Registers a new derived Parameter given by the arguments in the ScaleAndBlock File.
  //Currently only supports IFNONZERO[ ... ] and SLHA block parameters
 void SARAHManager::RegisterDerivedParameters(std::vector<std::pair<std::string,std::string>> derivedparameters)
    {
    for( auto it = derivedparameters.begin(); it != derivedparameters.end(); it++) //first get all the parameter definitions
		{	
        if(((*it).second).find("IFNONZERO")==std::string::npos)
            {
           LesHouchesAccordBlockEntryManager::RegisterNewParameter(LesHouchesAccordBlockEntryManager::CreateNewBlockEntry((*it).second),(*it).first);
            }
		}
    for( auto it = derivedparameters.begin(); it != derivedparameters.end(); it++) //TwoSource
      {
        if(((*it).second).find("IFNONZERO")!=std::string::npos)
            {
			std::string parnames(it->second);
			parnames = parnames.substr(parnames.find('[')+1, (parnames.find(']') - parnames.find('[') - 1));
			std::pair< bool, size_t > ifnonzero1 = LesHouchesAccordBlockEntryManager::RegisterParameter(parnames.substr(0,parnames.find(',')));
            std::pair< bool, size_t > ifnonzero2 = LesHouchesAccordBlockEntryManager::RegisterParameter(parnames.substr(parnames.find(',') + 1));
            if(ifnonzero1.first && ifnonzero2.first)
				{
                AddNewDerivedParameter((*it).first, new LhaTwoSourceFunctionoid(numberOfDistinctActiveParameters,ifnonzero1.second,ifnonzero2.second));

				}
            else
				{
                std::stringstream errorBuilder;
                errorBuilder
                << "RegisterDerivedParameter does not recognize the given IFNONZERO parameters, please check if the given parameters are correctly defined." << std::endl;
                throw std::runtime_error( errorBuilder.str() );
				}
            }
      }
	}
} /* namespace VevaciousPlusPlus */
