/*
 * GradientFromStartingPoints.cpp
 *
 *  Created on: Jun 30, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialMinimization/GradientFromStartingPoints.hpp"

namespace VevaciousPlusPlus
{

    GradientFromStartingPoints::GradientFromStartingPoints(
            PotentialFunction& potentialFunction,
            std::unique_ptr<StartingPointFinder> startingPointFinder,
            std::unique_ptr<GradientMinimizer> gradientMinimizer,
            double const extremumSeparationThresholdFraction,
            double const nonDsbRollingToDsbScalingFactor, bool global_Is_Panic ) :
            PotentialMinimizer( potentialFunction ),
            startingPointFinder( std::move(startingPointFinder) ),
            gradientMinimizer( std::move(gradientMinimizer) ),
            startingPoints(),
            extremumSeparationThresholdFraction( extremumSeparationThresholdFraction ),
            nonDsbRollingToDsbScalingFactor( nonDsbRollingToDsbScalingFactor ),
            global_Is_Panic(global_Is_Panic)
    {
        // This constructor is just an initialization list.
    }

    GradientFromStartingPoints::~GradientFromStartingPoints()
    {
    }


    // This first sets dsbVacuum from the input recorded in
    // potentialFunction.DsbFieldValues() using gradientMinimizer, then uses
    // startingPointFinder to find the starting points for gradientMinimizer,
    // then uses gradientMinimizer to minimize the potential at a temperature
    // given by minimizationTemperature, recording the found minima in
    // foundMinima. It also records the minima lower than dsbVacuum in
    // panicVacua, and of those, it sets panicVacuum to be the minimum in
    // panicVacua closest to dsbVacuum or the global minimum depending on what the user
    // set for global_Is_Panic. The default is the former.
    void GradientFromStartingPoints::FindMinima(
            double const minimizationTemperature )
    {
        gradientMinimizer->SetTemperature( minimizationTemperature );
        std::cout
                << std::endl
                << "DSB vacuum input: "
                << potentialFunction.FieldConfigurationAsMathematica(
                        potentialFunction.DsbFieldValues() );
        std::cout << std::endl;
        dsbVacuum = (*gradientMinimizer)( potentialFunction.DsbFieldValues() );
        std::cout
                << "Rolled to: "
                << dsbVacuum.AsMathematica( potentialFunction.FieldNames() );
        std::cout << std::endl;

        (*startingPointFinder)( startingPoints );
        double const
                thresholdSepartionSquared( ( extremumSeparationThresholdFraction
                                             * extremumSeparationThresholdFraction
                                             * dsbVacuum.LengthSquared() ) + 1.0 );
        double const thresholdSepartion( sqrt( thresholdSepartionSquared ) );
        PotentialMinimum foundMinimum;
        std::cout
                << std::endl
                << "Gradient-based minimization from a set of starting points:";
        for( std::vector< std::vector< double > >::const_iterator
                realSolution( startingPoints.begin() );
                realSolution < startingPoints.end();
        ++realSolution)
        {
            std::cout
                    << std::endl
                    << "Starting point: "
                    << potentialFunction.FieldConfigurationAsMathematica( *realSolution );
            std::cout << std::endl;
            foundMinimum = (*gradientMinimizer)( *realSolution );
            std::cout
                    << "Rolled to: "
                    << foundMinimum.AsMathematica( potentialFunction.FieldNames() );
            std::cout << std::endl;
            bool rolledToDsbOrSignFlip( ( foundMinimum.SquareDistanceTo( dsbVacuum )
                                          < thresholdSepartionSquared )
                                        ||
                                        !( IsNotPhaseRotationOfDsbVacuum( foundMinimum,
                                                                          thresholdSepartion ) ) );

            // We check to see if a starting point that was not the DSB minimum
            // rolled to the DSB minimum: if so, we scale the starting point's fields
            // by a factor and pass the scaled set of field values to
            // gradientMinimizer to roll, and then carry on based on this new
            // minimum. (We discovered in explorations with Vevacious 1 that
            // it could happen that the basin of attraction of the DSB minimum at
            // 1-loop level could grow so large that it would encompass tree-level
            // minima that belong in some sense to other 1-loop minima, which moved
            // very far away due to loop corrections, so even though their basins of
            // attraction also grew very large in the same way that of the DSB
            // minimum did, they moved enough that their tree-level minima were left
            // out.)
            if( rolledToDsbOrSignFlip
                &&
                ( dsbVacuum.SquareDistanceTo( *realSolution )
                  > thresholdSepartionSquared ) )
            {
                // We don't want to bother re-rolling the field origin, so we keep note
                // of how far away from the field origin realSolution is.
                double lengthSquared( 0.0 );
                std::vector< double > scaledPoint( *realSolution );
                for( std::vector< double >::iterator
                             scaledField( scaledPoint.begin() );
                     scaledField < scaledPoint.end();
                     ++scaledField )
                {
                    lengthSquared += ( (*scaledField) * (*scaledField) );
                    *scaledField *= nonDsbRollingToDsbScalingFactor;
                }
                if( lengthSquared > thresholdSepartionSquared )
                {
                    std::cout
                            << "Non-DSB-minimum starting point rolled to the DSB minimum, or a"
                            << " phase rotation, using the full potential. Trying a scaled"
                            << " starting point: "
                            << potentialFunction.FieldConfigurationAsMathematica( scaledPoint );
                    std::cout << std::endl;

                    foundMinimum = (*gradientMinimizer)( scaledPoint );
                    rolledToDsbOrSignFlip = ( foundMinimum.SquareDistanceTo( dsbVacuum )
                                              < thresholdSepartionSquared )
                                            ||
                                            !( IsNotPhaseRotationOfDsbVacuum( foundMinimum,
                                                                              thresholdSepartion ) );
                    std::cout
                            << "Rolled to: "
                            << foundMinimum.AsMathematica( potentialFunction.FieldNames() );
                    std::cout << std::endl;
                }
            }

            foundMinima.push_back( foundMinimum );

            // Here I do some checks so that we know minuit is behaving properly
            // This can be an issue with pathological parameter points or minima that
            // are too far away from the DSB, for which it makes no sense to calculate
            // tunneling without RG running anyway

            if(isnan(foundMinimum.FunctionValue()) || isnan(foundMinimum.FunctionError()) )
            {
                std::stringstream errorBuilder;
                errorBuilder << "Problem with Minuit, NaN given in minimum value/error. ";
                throw std::runtime_error( errorBuilder.str() );
            }

            // Here we set the panic vacuum to the closest lower minimum to the DSB minimum


            if(global_Is_Panic)
            {
                // Here we set the panic vacuum to the global minimum
                if( ( ( foundMinimum.FunctionValue() + foundMinimum.FunctionError() )
                      < dsbVacuum.FunctionValue() )
                    &&
                    !rolledToDsbOrSignFlip )
                {

                    panicVacuum = foundMinimum;

                    panicVacua.push_back( foundMinimum );
                }

            }else {

                if (((foundMinimum.FunctionValue() + foundMinimum.FunctionError())
                     < dsbVacuum.FunctionValue())
                    &&
                    !rolledToDsbOrSignFlip) {
                    if (panicVacua.empty()
                        ||
                        (foundMinimum.SquareDistanceTo(dsbVacuum)
                         < panicVacuum.SquareDistanceTo(dsbVacuum))) {
                        panicVacuum = foundMinimum;
                    }
                    panicVacua.push_back(foundMinimum);
                }

            }

        }

        std::cout
                << std::endl
                << "DSB vacuum = "
                << dsbVacuum.AsMathematica( potentialFunction.FieldNames() ) << std::endl;
        if( panicVacua.empty() )
        {
            std::cout
                    << "DSB vacuum is stable as far as the model file allows." << std::endl;
        }
        else
        {
            std::cout << "There are "
                      << panicVacua.size()
                      <<" panic vacua."
                      << std::endl;
            std::cout << "Panic vacuum used in tunneling = "
                      << panicVacuum.AsMathematica( potentialFunction.FieldNames() )
                      << std::endl;
        }
        std::cout << std::endl;
        std::cout << std::endl;
    }

} /* namespace VevaciousPlusPlus */
