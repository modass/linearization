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
 * Created by Stefan Schupp <stefan.schupp@tuwien.ac.at> on 23.11.22.
 */

#ifndef LINEARIZATION_SIMPLE_EXAMPLE_MATRICES_H
#define LINEARIZATION_SIMPLE_EXAMPLE_MATRICES_H

#include <hypro/types.h>

namespace linearization {

hypro::matrix_t<double> order3() {
	/* A:
	[[ 0. -1.  0.  0.  0.  0.]
	[ 1.  1.  0.  0.  0.  0.]
	[ 0. -1.  1. -2.  0.  3.]
	[ 0.  0.  2.  2.  0.  0.]
	[ 0.  0. -1.  0.  1.  0.]
	[ 0.  0. -1.  0.  0.  0.]]
	 */
	hypro::matrix_t<double> A = hypro::matrix_t<double>::Zero( 7, 7 );
	A( 0, 1 ) = -1.0;
	A( 1, 0 ) = 1.0;
	A( 1, 1 ) = 1.0;
	A( 2, 1 ) = -1.0;
	A( 2, 2 ) = 1.0;
	A( 2, 3 ) = -2.0;
	A( 2, 5 ) = 3.0;
	A( 3, 2 ) = 2.0;
	A( 3, 3 ) = 2.0;
	A( 4, 2 ) = -1.0;
	A( 4, 4 ) = 1.0;
	A( 5, 2 ) = -1.0;
	return A;
}

hypro::matrix_t<double> order5() {
	/* A:
	[[ 0. -1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
	 [ 1.  1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
	 [ 0. -1.  1. -2.  0.  3.  0.  0.  0.  0.  0.  0.  0.  0.]
	 [ 0.  0.  2.  2.  0.  0.  0.  0. -3.  0.  0.  0.  0.  0.]
	 [ 0.  0. -1.  0.  1.  0.  0. -2.  0.  5.  0.  0.  0.  0.]
	 [ 0.  0. -1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
	 [ 0.  0.  0.  0. -1.  0.  1.  0.  0.  0.  0. -2.  7.  0.]
	 [ 0.  0.  0. -2.  4.  0.  0.  2.  0.  0. -3.  0.  0.  0.]
	 [ 0.  0.  0.  1.  0.  0.  0.  0.  3.  0.  0.  0.  0.  0.]
	 [ 0.  0.  0.  0. -1.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
	 [ 0.  0.  0.  0.  0.  0.  0.  3. -3.  0.  3.  0.  0.  0.]
	 [ 0.  0.  0.  0.  0.  0.  6. -2.  0.  0.  0.  2.  0.  0.]
	 [ 0.  0.  0.  0.  0.  0. -1.  0.  0.  0.  0.  0.  0.  0.]
	 [ 0.  0.  0.  0.  0.  0. -1.  0.  0.  0.  0.  0.  0.  1.]]
	 */
	hypro::matrix_t<double> A = hypro::matrix_t<double>::Zero( 14, 14 );
	A << 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		  1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		  0.0, -1.0, 1.0, -2.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		  0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		  0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, -2.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0,
		  0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		  0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -2.0, 7.0, 0.0,
		  0.0, 0.0, 0.0, -2.0, 4.0, 0.0, 0.0, 2.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0,
		  0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		  0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, -3.0, 0.0, 3.0, 0.0, 0.0, 0.0,
		  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0, -2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0,
		  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
	return A;
}

}  // namespace linearization

#endif	// LINEARIZATION_SIMPLE_EXAMPLE_MATRICES_H
