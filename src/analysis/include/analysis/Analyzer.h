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

#ifndef LINEARIZATION_ANALYZER_H
#define LINEARIZATION_ANALYZER_H

#include "types.h"
#include <hypro/algorithms/reachability/Reach.h>
#include <hypro/datastructures/HybridAutomaton/HybridAutomaton.h>

namespace analysis
{

/**
 * Class that wraps reachability analysis for linear hybrid systems implemented
 * in HyPro
 */
template <typename Representation, typename Dynamics>
class Analyzer {
	using Number = typename Representation::NumberType;

  public:
	/// Constructor from an automaton and a set of initial valuations
	Analyzer( Dynamics&& dynamics, const std::vector<hypro::Condition<Number>>& initialValuations );
	/// Performs the analysis
	void run();

  protected:
	Dynamics mDynamics;
	std::vector<hypro::Condition<Number>> mInitialConditions;
	hypro::HybridAutomaton<Number> mAutomaton;
};

} // namespace analysis

#include "Analyzer.tpp"

#endif // LINEARIZATION_ANALYZER_H
