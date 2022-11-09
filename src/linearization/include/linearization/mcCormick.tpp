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
	// find best linearization points according to some heuristic (default: random) -> put this into the settings
	// lines 128-193
	// do this for lower and upper approximation
	// get: good linearization points for lower and upper approximation
	auto linearizationPoint = getLinearizationPoint( dynamics, settings );
	spdlog::trace( "Linearization point: {}", linearizationPoint );

	// do the relaxation in the chosen points

	// write result

	// create sample points
	auto combinator = hypro::Combinator( settings.subdivisions, settings.subdivisions.size() );
	while ( !combinator.end() ) {
		// returns a list of indices, which toll us which bucket of the discretization to use
		auto sample_indices = combinator();
		// create actual sample point and relaxation-objects
		std::vector<double> sample;
		std::vector<MC> relaxations;
		for ( std::size_t i = 0; i < sample_indices.size(); ++i ) {
			const auto& dom = settings.domain.intervals[i];
			sample.push_back( dom.l() + sample_indices[i] * ( mc::diam( dom ) / settings.subdivisions[i] ) );
			relaxations.push_back( MC( dom, sample.back() ) );
		}
		spdlog::trace( "Sample: {}", sample );
		MC relaxation = dynamics( relaxations );
	}

	throw utility::NotImplemented();
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
