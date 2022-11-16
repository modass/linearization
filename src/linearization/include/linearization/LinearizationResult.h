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
 * Created by Stefan Schupp <stefan.schupp@tuwien.ac.at> on 20.10.22.
 */

#ifndef LINEARIZATION_LINEARIZATIONRESULT_H
#define LINEARIZATION_LINEARIZATIONRESULT_H

#include <hypro/datastructures/HybridAutomaton/Condition.h>
#include <hypro/datastructures/HybridAutomaton/HybridAutomaton.h>

namespace linearization {

template <typename Number>
struct LinearizationResult {
	hypro::Condition<Number> initialCondition;	///< Polyhedral description of the linearized initial set
												// std::vector< dynamics;	///< Single-location automaton which holds the linearized dynamics
};

}  // namespace linearization

#endif	// LINEARIZATION_LINEARIZATIONRESULT_H
