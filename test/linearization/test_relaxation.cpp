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

#include <gtest/gtest.h>
#include <linearization/mcCormick.h>
#include <linearization/types.h>

namespace test::impl {

/*
template <typename N>
 struct constantMonomial {
	N operator()( const std::vector<N>& in ) const {
		N res = 0 * in[0] + 0 * in[1] + N( 2 );
		return res;
	}
};


template <typename N>
struct quadraticMonomial {
	// monomial: y = x^2, or: 0 = x^2 - y
	N operator()( const std::vector<N>& in ) const {
		N res = pow( in[0], 2 ) - in[1];
		return res;
	}
};
 */

}  // namespace test::impl

TEST( Relaxation, ConstantFunction ) {
	using namespace linearization;
	spdlog::set_level( spdlog::level::trace );

	std::function<MC( std::vector<MC> )> constantMonomial = []( const std::vector<MC>& in ) {MC res = 0 * in[0] + 0 * in[1] + MC( 2 );
		return res; };

	Domain d{ { Interval{ 0, 10 }, Interval{ 0, 10 } } };
	std::vector<std::size_t> subdivisions{ 5, 5 };
	Settings s{ d, subdivisions };
	LinearizationResult<double> result;
	// call linearization, the result should contain two constraints which under- and over-approximate the monomial within the domain
	linearizeMonomial( constantMonomial, s, result );
	// checks
	EXPECT_TRUE( result.initialCondition.size() == 1 );
	EXPECT_EQ( 2, result.initialCondition.getMatrix().rows() );
	EXPECT_EQ( 2, result.initialCondition.getMatrix().cols() );
	EXPECT_EQ( 2, result.initialCondition.getVector().rows() );
}

TEST( Relaxation, QuadraticFunction ) {
	using namespace linearization;
	spdlog::set_level( spdlog::level::debug );

	std::function<MC( std::vector<MC> )> quadraticMonomial = []( const std::vector<MC>& in ) {MC res = pow( in[0], 2 ) - in[1];
		return res; };

	Domain d{ { Interval{ -2, 2 }, Interval{ 0, 4 } } };
	std::vector<std::size_t> subdivisions{ 5, 5 };
	Settings s{ d, subdivisions };

	// test relaxation in specific points
	{
		auto res = getRelaxationInPoint( quadraticMonomial, d, Point{ 0, 0 } );
		EXPECT_EQ( 0, res.cv() );
		EXPECT_EQ( 4, res.cc() );
		EXPECT_EQ( 0, res.cvsub( 0 ) );
		EXPECT_EQ( -1, res.cvsub( 1 ) );
	}
	{
		auto res = getRelaxationInPoint( quadraticMonomial, d, Point{ 1, 0 } );
		EXPECT_EQ( 1, res.cv() );
		EXPECT_EQ( 4, res.cc() );
		EXPECT_EQ( 2, res.cvsub( 0 ) );
		EXPECT_EQ( -1, res.cvsub( 1 ) );
	}

	// test randomized relaxation
	LinearizationResult<double> result;
	// call linearization, the result should contain two constraints which under- and over-approximate the monomial within the domain
	linearizeMonomial( quadraticMonomial, s, result );
	// checks
	EXPECT_TRUE( result.initialCondition.size() == 1 );
	EXPECT_EQ( 2, result.initialCondition.getMatrix().rows() );
	EXPECT_EQ( 2, result.initialCondition.getMatrix().cols() );
	EXPECT_EQ( 2, result.initialCondition.getVector().rows() );
	// test some points
	int resolution = 100;
	for ( int i = 0; i <= resolution; ++i ) {
		double x = -2 + i * ( ( d.intervals[0].u() - d.intervals[0].l() ) / resolution );
		double y = pow( x, 2 );
		EXPECT_TRUE( result.initialCondition.contains( Point( { x, y } ) ) );
	}
}
