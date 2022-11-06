/*
 * Copyright (c) 2022 Stefan Schupp.
 * This file is part of the linearization project.
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/*
 * Created by Stefan Schupp <stefan.schupp@tuwien.ac.at> on 04.11.22.
 */

#include "linearization/heuristics/uniformRandom.h"

#include <random>

namespace linearization::impl {
Point generateRandomPoint( const Domain& domain ) {
	Point res = Point( hypro::vector_t<double>( domain.intervals.size() ) );
	std::random_device rd;
	std::mt19937 gen( rd() );
	Index pos{ 0 };
	for ( const auto& interval : domain.intervals ) {
		std::uniform_real_distribution<> dist( interval.l(), interval.u() );
		res[pos++] = dist( gen );
	}
	return res;
}
}  // namespace linearization::impl
