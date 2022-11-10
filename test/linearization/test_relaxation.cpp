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
 * Created by Stefan Schupp <stefan.schupp@tuwien.ac.at> on 04.11.22.
 */

#include <gtest/gtest.h>
#include <linearization/mcCormick.h>
#include <linearization/types.h>

namespace test::impl {

template <typename N>
struct constantMonomial {
	N operator()( const std::vector<N>& in ) const {
		N res = 0 * in[0] + 0 * in[1] + N( 2 );
		return res;
	}
};

template <typename N>
struct quadraticMonomial {
	N operator()( const std::vector<N>& in ) const {
		N res = pow( in[0], 2 ) + 0 * in[1];
		return res;
	}
};

}  // namespace test::impl

TEST( Relaxation, ConstantFunction ) {
	using namespace linearization;
	spdlog::set_level( spdlog::level::trace );

	Domain d{ { Interval{ 0, 10 }, Interval{ 0, 10 } } };
	std::vector<std::size_t> subdivisions{ 5, 5 };
	Settings s{ d, subdivisions };
	LinearizationResult<double> result;
	// call linearization, the result should contain two constraints which under- and over-approximate the monomial within the domain
	linearizeMonomial( test::impl::constantMonomial<MC>(), s, result );
	// checks
	EXPECT_TRUE( result.initialCondition.size() == 1 );
	EXPECT_EQ( 2, result.initialCondition.getMatrix().rows() );
	EXPECT_EQ( 2, result.initialCondition.getMatrix().cols() );
	EXPECT_EQ( 2, result.initialCondition.getVector().rows() );
}

TEST( Relaxation, QuadraticFunction ) {
	using namespace linearization;
	spdlog::set_level( spdlog::level::debug );

	Domain d{ { Interval{ 0, 2 }, Interval{ 0, 2 } } };
	std::vector<std::size_t> subdivisions{ 5, 5 };
	Settings s{ d, subdivisions };
	LinearizationResult<double> result;
	// call linearization, the result should contain two constraints which under- and over-approximate the monomial within the domain
	linearizeMonomial( test::impl::quadraticMonomial<MC>(), s, result );
	// checks
	EXPECT_TRUE( result.initialCondition.size() == 1 );
	EXPECT_EQ( 2, result.initialCondition.getMatrix().rows() );
	EXPECT_EQ( 2, result.initialCondition.getMatrix().cols() );
	EXPECT_EQ( 2, result.initialCondition.getVector().rows() );
}
