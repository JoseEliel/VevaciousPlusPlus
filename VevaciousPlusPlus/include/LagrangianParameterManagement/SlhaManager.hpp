/*
 * SlhaManager.hpp
 *
 *  Created on: Apr 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHAMANAGER_HPP_
#define SLHAMANAGER_HPP_

#include "CommonIncludes.hpp"

namespace VevaciousPlusPlus
{
  // This class is an abstract base class to try to impose some consistency on
  // PotentialMinimizers and TunnelingCalculators so that they share the data
  // from an SLHA file.
  class SlhaManager : public LagrangianParameterManager
  {
  public:
    SlhaManager();
    virtual ~SlhaManager();


    void UpdateSlhaData( std::string const& slhaFilename );

    // This should return the lowest Q found among the SLHA blocks, returning
    // 0.0 if no block had a scale greater than 0.0 GeV.
    virtual double LowestBlockScale() const = 0;

    // This should return the highest Q found among the SLHA blocks, returning
    // 0.0 if no block had a scale greater than 0.0 GeV.
    virtual double HighestBlockScale() const = 0;


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
