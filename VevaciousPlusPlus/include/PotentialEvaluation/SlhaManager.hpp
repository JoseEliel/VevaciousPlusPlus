/*
 * SlhaManager.hpp
 *
 *  Created on: Apr 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHAMANAGER_HPP_
#define SLHAMANAGER_HPP_

#include "../StandardIncludes.hpp"
#include "BOLlib/include/BasicObserved.hpp"

namespace VevaciousPlusPlus
{
  // This class is an abstract base class to try to impose some consistency on
  // PotentialMinimizers and TunnelingCalculators so that they share the data
  // from an SLHA file.
  class SlhaManager : public BOL::BasicObserved
  {
  public:
    SlhaManager();
    virtual
    ~SlhaManager();


    void UpdateSlhaData( std::string const& slhaFilename );


  protected:
    virtual void ReadFile( std::string const& slhaFilename ) = 0;
  };




  inline void SlhaManager::UpdateSlhaData( std::string const& slhaFilename )
  {
    ReadFile( slhaFilename );
    updateObservers();
  }
} /* namespace VevaciousPlusPlus */
#endif /* SLHAMANAGER_HPP_ */
