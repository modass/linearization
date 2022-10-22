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

#include "mcCormick.h"

#include <MCpp/include/interval.hpp>
#include <MCpp/include/mccormick.hpp>
#include <hypro/util/sequenceGeneration/SequenceGenerator.h>
#include <utility/Exceptions.h>

namespace linearization {

template <typename Function>
LinearizationResult<double> linearize( const Function& dynamics, const Settings& settings ) {
	using Interval = mc::Interval;
	using MC = mc::McCormick<Interval>;

	LinearizationResult<double> res;

	MC::options.ENVEL_USE = true;
	MC::options.ENVEL_MAXIT = 100;
	MC::options.ENVEL_TOL = 1e-12;
	MC::options.MVCOMP_USE = true;

	// create sample points
	auto combinator = hypro::Combinator( settings.subdivisions, settings.subdivisions.size() );
	while ( !combinator.end() ) {
		// returns a list of indices, which toll us which bucket of the discretization to use
		auto sample_indices = combinator();
		// create actual sample point and relaxation-objects
		std::vector<double> sample;
		std::vector<MC> relaxations;
		for ( std::size_t i = 0; i < sample_indices.size(); ++i ) {
			const auto& dom = settings.domain[i];
			sample.push_back( dom.l() + sample_indices[i] * ( diam( dom ) / settings.subdivisions[i] ) );
			relaxations.push_back( MC( dom, sample.back() ) );
		}
		MC relaxation = dynamics( relaxations );
	}

	throw utility::NotImplemented();
	return res;
}

}  // namespace linearization