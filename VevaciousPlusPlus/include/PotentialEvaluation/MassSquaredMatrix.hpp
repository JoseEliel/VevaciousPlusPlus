/*
 * MassSquaredMatrix.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MASSSQUAREDMATRIX_HPP_
#define MASSSQUAREDMATRIX_HPP_

#include "../StandardIncludes.hpp"
#include "Eigen/Dense"

namespace VevaciousPlusPlus
{
  class MassSquaredMatrix
  {
  public:
    MassSquaredMatrix();
    virtual
    ~MassSquaredMatrix();



  protected:
    Eigen::MatrixXd eigenMatrix;
  };

} /* namespace VevaciousPlusPlus */
#endif /* MASSSQUAREDMATRIX_HPP_ */
