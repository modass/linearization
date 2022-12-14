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

#ifndef LINEARIZATION_TYPES_H
#define LINEARIZATION_TYPES_H

#include <MCpp/include/interval.hpp>
#include <MCpp/include/mccormick.hpp>
#include <hypro/datastructures/Point.h>
#include <vector>

namespace linearization {

using Point = hypro::Point<double>;
using Index = Eigen::Index;

using Interval = mc::Interval;
using MC = mc::McCormick<Interval>;

template <typename N, typename... Rargs>
using Function = std::function<N( Rargs... )>();

struct Domain {
	std::vector<Interval> intervals;
};

enum class Approximation {
	OVER = 0,
	UNDER = 1
};

}  // namespace linearization

#endif	// LINEARIZATION_TYPES_H
