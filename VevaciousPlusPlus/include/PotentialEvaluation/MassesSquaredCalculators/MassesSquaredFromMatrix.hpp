/*
 * MassesSquaredFromMatrix.hpp
 *
 *  Created on: Oct 09, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MASSESSQUAREDFROMMATRIX_HPP_
#define MASSESSQUAREDFROMMATRIX_HPP_

#include "PotentialEvaluation/MassesSquaredCalculator.hpp"
#include "Eigen/Dense"
#include <map>
#include <string>
#include <vector>

namespace VevaciousPlusPlus
{
  template< typename ElementType > class MassesSquaredFromMatrix :
                                                 public MassesSquaredCalculator
  {
  public:
    typedef typename
    Eigen::Matrix< ElementType, Eigen::Dynamic, Eigen::Dynamic > EigenMatrix;

    MassesSquaredFromMatrix( size_t numberOfRows,
                    std::map< std::string, std::string > const& attributeMap );
    MassesSquaredFromMatrix( MassesSquaredFromMatrix const& copySource );
    MassesSquaredFromMatrix();
    virtual ~MassesSquaredFromMatrix();


    // This returns the eigenvalues of the matrix, using the values for the
    // Lagrangian parameters found in parameterValues and the values for the
    // fields found in fieldConfiguration.
    virtual std::vector< double >
    MassesSquared( std::vector< double > const& parameterValues,
                   std::vector< double > const& fieldConfiguration ) const;

    // This returns the eigenvalues of the matrix, using the values for the
    // Lagrangian parameters from the last call of UpdateForFixedScale and the
    // values for the fields found in fieldConfiguration.
    virtual std::vector< double >
    MassesSquared( std::vector< double > const& fieldConfiguration ) const;

    size_t NumberOfRows() const{ return numberOfRows; }


  protected:
    size_t numberOfRows;

    // This should return a matrix of the values of the elements for a field
    // configuration given by fieldConfiguration, using the values for the
    // Lagrangian parameters found in parameterValues.
    virtual EigenMatrix
    CurrentValues( std::vector< double > const& parameterValues,
                   std::vector< double > const& fieldConfiguration ) const = 0;

    // This should return a matrix of the values of the elements for a field
    // configuration given by fieldConfiguration, using the values for the
    // Lagrangian parameters from the last call of UpdateForFixedScale.
    virtual EigenMatrix
    CurrentValues( std::vector< double > const& fieldConfiguration ) const = 0;
  };



  template< typename ElementType > inline
  MassesSquaredFromMatrix< ElementType >::MassesSquaredFromMatrix(
                                                           size_t numberOfRows,
                   std::map< std::string, std::string > const& attributeMap ) :
    MassesSquaredCalculator( attributeMap ),
    numberOfRows( numberOfRows )
  {
    // This constructor is just an initialization list.
  }

  template< typename ElementType > inline
  MassesSquaredFromMatrix< ElementType >::MassesSquaredFromMatrix(
                   MassesSquaredFromMatrix< ElementType > const& copySource ) :
    MassesSquaredCalculator( copySource ),
    numberOfRows( copySource.numberOfRows )
  {
    // This constructor is just an initialization list.
  }

  template< typename ElementType > inline
  MassesSquaredFromMatrix< ElementType >::MassesSquaredFromMatrix() :
    MassesSquaredCalculator(),
    numberOfRows( 0 )
  {
    // This constructor is just an initialization list.
  }

  template< typename ElementType > inline
  MassesSquaredFromMatrix< ElementType >::~MassesSquaredFromMatrix()
  {
    // This does nothing.
  }

  // This returns the eigenvalues of the matrix, using the values for the
  // Lagrangian parameters found in parameterValues and the values for the
  // fields found in fieldConfiguration.
  template< typename ElementType > inline std::vector< double >
  MassesSquaredFromMatrix< ElementType >::MassesSquared(
                                  std::vector< double > const& parameterValues,
                        std::vector< double > const& fieldConfiguration ) const
  {
    Eigen::SelfAdjointEigenSolver< EigenMatrix >
    eigenvalueFinder( CurrentValues( parameterValues,
                                     fieldConfiguration ),
                      Eigen::EigenvaluesOnly );
    std::vector< double > massesSquared( numberOfRows );
    for( size_t rowIndex( 0 );
         rowIndex < numberOfRows;
         ++rowIndex )
    {
      massesSquared[ rowIndex ] = eigenvalueFinder.eigenvalues()( rowIndex );
    }
    return massesSquared;
  }

  // This returns the eigenvalues of the matrix, using the values for the
  // Lagrangian parameters from the last call of UpdateForFixedScale and the
  // values for the fields found in fieldConfiguration.
  template< typename ElementType > inline std::vector< double >
  MassesSquaredFromMatrix< ElementType >::MassesSquared(
                        std::vector< double > const& fieldConfiguration ) const
  {
    Eigen::SelfAdjointEigenSolver< EigenMatrix >
    eigenvalueFinder( CurrentValues( fieldConfiguration ),
                      Eigen::EigenvaluesOnly );
    std::vector< double > massesSquared( numberOfRows );
    for( size_t rowIndex( 0 );
         rowIndex < numberOfRows;
         ++rowIndex )
    {
      massesSquared[ rowIndex ] = eigenvalueFinder.eigenvalues()( rowIndex );
    }
    return massesSquared;
  }

} /* namespace VevaciousPlusPlus */
#endif /* MASSESSQUAREDFROMMATRIX_HPP_ */
