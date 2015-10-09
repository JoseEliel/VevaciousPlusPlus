/*
 * MassesSquaredFromMatrix.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef OLDMASSESSQUAREDFROMMATRIX_HPP_
#define OLDMASSESSQUAREDFROMMATRIX_HPP_

#include "CommonIncludes.hpp"
#include "../MassesSquaredCalculator.hpp"
#include "Eigen/Dense"

namespace VevaciousPlusPlus
{
  template< typename ElementType > class OldMassesSquaredFromMatrix :
                                              public OldMassesSquaredCalculator
  {
  public:
    typedef typename
    Eigen::Matrix< ElementType, Eigen::Dynamic, Eigen::Dynamic > EigenMatrix;

    OldMassesSquaredFromMatrix( size_t numberOfRows,
                    std::map< std::string, std::string > const& attributeMap );
    OldMassesSquaredFromMatrix( OldMassesSquaredFromMatrix const& copySource );
    OldMassesSquaredFromMatrix();
    virtual ~OldMassesSquaredFromMatrix();


    // This returns the eigenvalues of the matrix, with all functionoids
    // evaluated at the last scale which was used to update them.
    virtual std::vector< double >
    MassesSquared( std::vector< double > const& fieldConfiguration ) const;

    // This returns the eigenvalues of the matrix, with all functionoids
    // evaluated at the natural exponent of logarithmOfScale.
    virtual std::vector< double >
    MassesSquared( std::vector< double > const& fieldConfiguration,
                   double const logarithmOfScale ) const;

    size_t NumberOfRows() const{ return numberOfRows; }


  protected:
    size_t numberOfRows;

    // This should return a matrix of the values of the elements for a field
    // configuration given by fieldConfiguration, with all functionoids
    // evaluated at the last scale which was used to update them.
    virtual EigenMatrix
    CurrentValues( std::vector< double > const& fieldConfiguration ) const = 0;

    // This should return a matrix of the values of the elements for a field
    // configuration given by fieldConfiguration, with all functionoids
    // evaluated at the natural exponent of logarithmOfScale.
    virtual EigenMatrix
    CurrentValues( std::vector< double > const& fieldConfiguration,
                   double const logarithmOfScale ) const = 0;
  };



  template< typename ElementType > inline
  OldMassesSquaredFromMatrix< ElementType >::OldMassesSquaredFromMatrix(
                                                           size_t numberOfRows,
                   std::map< std::string, std::string > const& attributeMap ) :
    OldMassesSquaredCalculator( attributeMap ),
    numberOfRows( numberOfRows )
  {
    // This constructor is just an initialization list.
  }

  template< typename ElementType > inline
  OldMassesSquaredFromMatrix< ElementType >::OldMassesSquaredFromMatrix(
                   OldMassesSquaredFromMatrix< ElementType > const& copySource ) :
    OldMassesSquaredCalculator( copySource ),
    numberOfRows( copySource.numberOfRows )
  {
    // This constructor is just an initialization list.
  }

  template< typename ElementType > inline
  OldMassesSquaredFromMatrix< ElementType >::OldMassesSquaredFromMatrix() :
    OldMassesSquaredCalculator(),
    numberOfRows( 0 )
  {
    // This constructor is just an initialization list.
  }

  template< typename ElementType > inline
  OldMassesSquaredFromMatrix< ElementType >::~OldMassesSquaredFromMatrix()
  {
    // This does nothing.
  }

  // This returns the eigenvalues of the matrix, with all functionoids
  // evaluated at the last scale which was used to update them.
  template< typename ElementType > inline std::vector< double >
  OldMassesSquaredFromMatrix< ElementType >::MassesSquared(
                        std::vector< double > const& fieldConfiguration ) const
  {
    Eigen::SelfAdjointEigenSolver< EigenMatrix >
    eigenvalueFinder( CurrentValues( fieldConfiguration ),
                      Eigen::EigenvaluesOnly );
    std::vector< double > massesSquared( numberOfRows );
    for( size_t whichIndex( 0 );
         whichIndex < numberOfRows;
         ++whichIndex )
    {
      massesSquared[ whichIndex ]
      = eigenvalueFinder.eigenvalues()( whichIndex );
    }
    return massesSquared;
  }

  // This returns the eigenvalues of the matrix, with all functionoids
  // evaluated at the natural exponent of logarithmOfScale.
  template< typename ElementType > inline std::vector< double >
  OldMassesSquaredFromMatrix< ElementType >::MassesSquared(
                               std::vector< double > const& fieldConfiguration,
                                          double const logarithmOfScale ) const
  {
    Eigen::SelfAdjointEigenSolver< EigenMatrix >
    eigenvalueFinder( CurrentValues( fieldConfiguration,
                                     logarithmOfScale ),
                      Eigen::EigenvaluesOnly );
    std::vector< double > massesSquared( numberOfRows );
    for( size_t whichIndex( 0 );
         whichIndex < numberOfRows;
         ++whichIndex )
    {
      massesSquared[ whichIndex ]
      = eigenvalueFinder.eigenvalues()( whichIndex );
    }
    return massesSquared;
  }

} /* namespace VevaciousPlusPlus */
#endif /* OLDMASSESSQUAREDFROMMATRIX_HPP_ */
