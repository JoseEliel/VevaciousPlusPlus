/*
 * ComplexMassSquaredMatrix.cpp
 *
 *  Created on: Oct 09, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialEvaluation/MassesSquaredCalculators/ComplexMassSquaredMatrix.hpp"

namespace VevaciousPlusPlus
{

  ComplexMassSquaredMatrix::ComplexMassSquaredMatrix(
                                                     size_t const numberOfRows,
                   std::map< std::string, std::string > const& attributeMap ) :
    BaseComplexMassMatrix( numberOfRows,
                           attributeMap )
  {
    // This constructor is just an initialization list.
  }

  ComplexMassSquaredMatrix::ComplexMassSquaredMatrix(
                                 ComplexMassSquaredMatrix const& copySource ) :
    BaseComplexMassMatrix( copySource )
  {
    // This constructor is just an initialization list.
  }

  ComplexMassSquaredMatrix::ComplexMassSquaredMatrix() :
    BaseComplexMassMatrix()
  {
    // This constructor is just an initialization list.
  }

  ComplexMassSquaredMatrix::~ComplexMassSquaredMatrix()
  {
    // This does nothing.
  }

  // This returns a matrix of the values of the elements for a field
   // configuration given by fieldConfiguration, using the values for the
   // Lagrangian parameters found in parameterValues.
  Eigen::MatrixXcd ComplexMassSquaredMatrix::CurrentValues(
                                  std::vector< double > const& parameterValues,
                        std::vector< double > const& fieldConfiguration ) const
  {
    size_t rowsTimesLength( 0 );
    Eigen::MatrixXcd valuesMatrix( numberOfRows,
                                   numberOfRows );
    for( size_t rowIndex( 0 );
         rowIndex < numberOfRows;
         ++rowIndex )
    {
      for( size_t columnIndex( 0 );
           columnIndex < rowIndex;
           ++columnIndex )
      {
        valuesMatrix.coeffRef( rowIndex,
                               columnIndex ).real()
        = matrixElements[ rowsTimesLength + columnIndex ].first(
                                                               parameterValues,
                                                          fieldConfiguration );
        valuesMatrix.coeffRef( rowIndex,
                               columnIndex ).imag()
        = matrixElements[ rowsTimesLength + columnIndex ].second(
                                                               parameterValues,
                                                          fieldConfiguration );
        // The Eigen routines don't bother looking at elements of valuesMatrix
        // where columnIndex > rowIndex, so we don't even bother filling them
        // with the conjugates of the transpose.
      }
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex ).real()
      = matrixElements[ rowsTimesLength + rowIndex ].first( parameterValues,
                                                          fieldConfiguration );
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex ).imag() = 0.0;
      rowsTimesLength += numberOfRows;
    }
    return valuesMatrix;
  }

  // This returns a matrix of the values of the elements for a field
  // configuration given by fieldConfiguration, using the values for the
  // Lagrangian parameters from the last call of UpdateForFixedScale.
  Eigen::MatrixXcd ComplexMassSquaredMatrix::CurrentValues(
                        std::vector< double > const& fieldConfiguration ) const
  {
    size_t rowsTimesLength( 0 );
    Eigen::MatrixXcd valuesMatrix( numberOfRows,
                                   numberOfRows );
    for( size_t rowIndex( 0 );
         rowIndex < numberOfRows;
         ++rowIndex )
    {
      for( size_t columnIndex( 0 );
           columnIndex < rowIndex;
           ++columnIndex )
      {
        valuesMatrix.coeffRef( rowIndex,
                               columnIndex ).real()
        = matrixElements[ rowsTimesLength + columnIndex ].first(
                                                          fieldConfiguration );
        valuesMatrix.coeffRef( rowIndex,
                               columnIndex ).imag()
        = matrixElements[ rowsTimesLength + columnIndex ].second(
                                                          fieldConfiguration );
        // The Eigen routines don't bother looking at elements of valuesMatrix
        // where columnIndex > rowIndex, so we don't even bother filling them
        // with the conjugates of the transpose.
      }
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex ).real()
      = matrixElements[ rowsTimesLength + rowIndex ].first(
                                                          fieldConfiguration );
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex ).imag() = 0.0;
      rowsTimesLength += numberOfRows;
    }
    return valuesMatrix;
  }

} /* namespace VevaciousPlusPlus */
