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
struct monomial {
	N operator()( const std::vector<N>& in ) const {
		return N( 2 );
	}
};

}  // namespace test::impl

TEST( Relaxation, ConstantFunction ) {
	using namespace linearization;
	spdlog::set_level(spdlog::level::trace);

	Domain d{ { Interval{ 0, 10 }, Interval{ 0, 10 } } };
	std::vector<std::size_t> subdivisions{ 5, 5 };
	Settings s{ d, subdivisions };
	LinearizationResult<double> result;

	linearizeMonomial( test::impl::monomial<MC>(), s, result );
}
