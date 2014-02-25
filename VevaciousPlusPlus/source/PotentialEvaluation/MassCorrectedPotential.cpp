/*
 * MassCorrectedPotential.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  MassCorrectedPotential::MassCorrectedPotential(
                                           std::string const& modelFilename ) :
    HomotopyContinuationReadyPotential()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MassCorrectedPotential::MassCorrectedPotential( \""
    << modelFilename << "\" )";
    std::cout << std::endl;/**/
  }

  MassCorrectedPotential::~MassCorrectedPotential()
  {
    // This does nothing.
  }


  MassCorrectedPotential::MassCorrectedPotential()
  {
    // This protected constructor is just an initialization list only used by
    // derived classes which are going to fill up the data members in their own
    // constructors.
  }

} /* namespace VevaciousPlusPlus */
