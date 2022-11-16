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
 * Created by Stefan Schupp <stefan.schupp@tuwien.ac.at> on 03.11.22.
 */

#include <gtest/gtest.h>
#include <linearization/heuristics/uniformRandom.h>

TEST( Heuristics, UniformRandom ) {
	linearization::Domain dom{ { linearization::Interval{ 0, 10 }, linearization::Interval{ 0, 10 }, linearization::Interval{ 0, 10 } } };

	auto contains = []( const linearization::Domain& d, const linearization::Point& p ) {
		for ( int k = 0; k < d.intervals.size(); ++k ) {
			if ( p[k] < d.intervals[k].l() || p[k] > d.intervals[k].u() ) {
				return false;
			}
		}
		return true;
	};

	// fuzzy testing
	for ( int k = 0; k < 100; ++k ) {
		auto sample = linearization::impl::generateRandomPoint( dom );
		EXPECT_TRUE( contains( dom, sample ) );
	}

	SUCCEED();
}
