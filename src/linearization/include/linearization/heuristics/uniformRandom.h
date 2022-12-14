/*
 * Copyright (c) 2022 Stefan Schupp.
 * This file is part of the linearization project.
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
 *
 */

/*
 * Created by Stefan Schupp <stefan.schupp@tuwien.ac.at> on 03.11.22.
 */

#ifndef LINEARIZATION_UNIFORMRANDOM_H
#define LINEARIZATION_UNIFORMRANDOM_H

#include "../mcCormick.h"
#include "../types.h"

#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>

namespace linearization {

namespace impl {
Point generateRandomPoint( const Domain& domain );
}  // namespace impl

template <typename Function>
Point getLinearizationPoint( const Function& dynamics, const Settings& settings, Approximation mode ) {
	assert( settings.subdivisions.size() == settings.domain.intervals.size() );
	// solution candidates
	double best_error = std::numeric_limits<double>::infinity();
	Point res;

	// make sure mcpp-settings are set
	settings.setMcppSettings();

	const auto& domain = settings.domain;
	for ( int iteration = 0; iteration < 10; ++iteration ) {
		auto dim = domain.intervals.size();
		/*set the point for which to compute the affine approximation (linearization point=(x_ref,y_ref)), sampling from a uniform distribution. The support for the considered uniform distribution is equal
		to the initial interval given for each coordinate.*/
		auto sample = impl::generateRandomPoint( domain );

		// McCormick relaxations of function myfunct on [XL,XU]\times[YL,YU] at (Xrel, Yrel):
		MC sample_relaxation = getRelaxationInPoint( dynamics, domain, sample );

		double squared_error = 0;
		// compute relaxation for grid
		auto combinator = hypro::Combinator( settings.subdivisions, settings.subdivisions.size() );
		while ( !combinator.end() ) {
			// returns a list of indices, which toll us which bucket of the discretization to use
			auto gridpoint_indices = combinator();
			// create actual grid point
			auto gridpoint = Point( hypro::vector_t<double>( dim ) );
			for ( std::size_t j = 0; j < gridpoint_indices.size(); ++j ) {
				const auto& dom = settings.domain.intervals[j];
				gridpoint[j] = ( dom.l() + gridpoint_indices[j] * ( mc::diam( dom ) / settings.subdivisions[j] ) );
			}
			MC gridpoint_relaxation = getRelaxationInPoint( dynamics, domain, gridpoint );
			// add error
			double local_error = 0;
			if ( mode == linearization::Approximation::UNDER ) {
				for ( std::size_t j = 0; j < gridpoint_indices.size(); ++j ) {
					local_error += sample_relaxation.cvsub( j ) * ( gridpoint[j] - sample[j] );
				}
				squared_error += pow( gridpoint_relaxation.cv() - sample_relaxation.cv() + local_error, 2 );
			} else {
				for ( std::size_t j = 0; j < gridpoint_indices.size(); ++j ) {
					local_error += sample_relaxation.ccsub( j ) * ( gridpoint[j] - sample[j] );
				}
				squared_error += pow( gridpoint_relaxation.cc() - sample_relaxation.cv() + local_error, 2 );
			}
		}
		double error = std::sqrt( squared_error );

		/*if the MSE obtained generating the affine approximation at the current linearization point leads to a lower MSE than the minimum obtained so far, then updates the final linearization point and the minimum MSE.*/
		if ( error < best_error ) {
			res = sample;
			best_error = error;
		}
	}
	if ( best_error == std::numeric_limits<double>::infinity() ) {
		throw std::logic_error( "Could not find a suitable point for relaxation" );
	}
	return res;
}
}  // namespace linearization

#endif	// LINEARIZATION_UNIFORMRANDOM_H
