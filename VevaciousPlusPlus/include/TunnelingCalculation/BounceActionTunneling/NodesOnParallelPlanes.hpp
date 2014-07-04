/*
 * NodesOnParallelPlanes.hpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef NODESONPARALLELPLANES_HPP_
#define NODESONPARALLELPLANES_HPP_

#include "CommonIncludes.hpp"
#include "NodesOnPlanes.hpp"

namespace VevaciousPlusPlus
{

  class NodesOnParallelPlanes
  {
  public:
    NodesOnParallelPlanes();
    virtual ~NodesOnParallelPlanes();


    // This sets up all the nodes based on the numbers given in
    // pathParameterization as a sequence of nodes, repeatedly calling
    // SetOneNode with a subvector of pathParameterization.
    virtual void
    SetAllNodes( std::vector< double > const& pathParameterization );

  protected:
  };

} /* namespace VevaciousPlusPlus */
#endif /* NODESONPARALLELPLANES_HPP_ */
