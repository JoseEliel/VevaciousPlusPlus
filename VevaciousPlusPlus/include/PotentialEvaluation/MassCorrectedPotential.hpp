/*
 * MassCorrectedPotential.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MASSCORRECTEDPOTENTIAL_HPP_
#define MASSCORRECTEDPOTENTIAL_HPP_

#include "../StandardIncludes.hpp"
#include "HomotopyContinuationReadyPotential.hpp"

namespace VevaciousPlusPlus
{
  class MassCorrectedPotential : public HomotopyContinuationReadyPotential
  {
  public:
    MassCorrectedPotential( std::string const& modelFilename );
    virtual
    ~MassCorrectedPotential();


  protected:
    MassCorrectedPotential();
  };

} /* namespace VevaciousPlusPlus */
#endif /* MASSCORRECTEDPOTENTIAL_HPP_ */
