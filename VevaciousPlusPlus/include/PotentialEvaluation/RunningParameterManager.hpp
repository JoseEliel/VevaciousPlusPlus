/*
 * RunningParameterManager.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef RUNNINGPARAMETERMANAGER_HPP_
#define RUNNINGPARAMETERMANAGER_HPP_

#include "../StandardIncludes.hpp"
#include "ParameterFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class RunningParameterManager
  {
  public:
    RunningParameterManager();
    virtual
    ~RunningParameterManager();


    // This adds to the list of SLHA blocks that will be read in when updating
    // running parameters. An alias for the block can be provided.
    void AddValidSlhaBlock( std::string const& slhaBlock,
                            std::string const& aliasString );

    // This creates a new ParameterFunctionoid and a string mapping to it.
    void CreateDerivedParameter( std::string const& parameterString );

    // This returns a pointer to the ParameterFunctionoid mapped to by
    // parameterString, returning NULL if there is none.
    ParameterFunctionoid* GetFunctionoid( std::string const& parameterString );


  protected:
  };

} /* namespace VevaciousPlusPlus */
#endif /* RUNNINGPARAMETERMANAGER_HPP_ */
