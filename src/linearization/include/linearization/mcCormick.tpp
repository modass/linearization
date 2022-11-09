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

//
// Created by Stefan Schupp on 17.10.22.
//

#pragma once

#include "heuristics/uniformRandom.h"
#include "mcCormick.h"
#include "types.h"

#include <MCpp/include/mccormick.hpp>
#include <hypro/util/sequenceGeneration/SequenceGenerator.h>
#include <spdlog/fmt/bundled/ostream.h>
#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/spdlog.h>
#include <utility/Exceptions.h>

namespace linearization {

template <typename Function>
void linearizeMonomial( const Function& dynamics, const Settings& settings, LinearizationResult<double>& result ) {
	auto dim = settings.domain.intervals.size();
	// find best linearization points according to some heuristic (default: random) -> put this into the settings
	// lines 128-193
	// do this for lower and upper approximation
	// get: good linearization points for lower and upper approximation
	auto linearizationPoint = getLinearizationPoint( dynamics, settings );
	spdlog::trace( "Linearization point: {}", linearizationPoint );

	// do the relaxation in the chosen points
	auto rel = getRelaxationInPoint( dynamics, settings.domain, linearizationPoint );

	// write result
	hypro::matrix_t<double> A = hypro::matrix_t<double>::Zero( 2, dim );
	hypro::vector_t<double> b = hypro::vector_t<double>::Zero( 2 );

	b( 0 ) = -rel.cv();
	b( 1 ) = -rel.cc();
	for ( int k = 0; k < dim; ++k ) {
		A( 0, k ) = rel.cvsub( k );
		A( 1, k ) = rel.ccsub( k );
		b( 0 ) += rel.cvsub( k ) * linearizationPoint[k];
		b( 1 ) += rel.ccsub( k ) * linearizationPoint[k];
	}
	result.initialCondition = hypro::Condition<double>( A, b );
}

template <typename Function>
MC getRelaxationInPoint( const Function& dynamics, const Domain& domain, const Point& point ) {
	auto dim = domain.intervals.size();
	std::vector<MC> relaxations;
	for ( std::size_t i = 0; i < dim; ++i ) {
		// create relaxation in point per dimension
		relaxations.push_back( MC( domain.intervals[i], point[i] ) );
		// defining subgradient components
		relaxations.back().sub( dim, i );
	}
	spdlog::trace( "Relaxations pre-execution of the function: {}", relaxations );
	// compute the McCormick convex and concave relaxations of myfunc at (Xrel,Yrel) along with subgradients of these relaxations.
	MC res = dynamics( relaxations );
	spdlog::debug( "Created relaxation {} in point {}", res, point );
	return res;
}

}  // namespace linearization
