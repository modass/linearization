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
 * Created by Stefan Schupp <stefan.schupp@tuwien.ac.at> on 21.11.22.
 */

#ifndef LINEARIZATION_CONVERSION_H
#define LINEARIZATION_CONVERSION_H

#include "types.h"

namespace linearization {

template <typename Number>
inline carl::Interval<double> convert( const Interval& in ) {
	return carl::Interval<double>( in.l(), in.u() );
}

}  // namespace linearization

#endif	// LINEARIZATION_CONVERSION_H
