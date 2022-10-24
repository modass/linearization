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

#include <analysis/Analyzer.h>
#include <gtest/gtest.h>
#include <hypro/representations/GeometricObjectBase.h>

// Demonstrate some basic assertions.
TEST( AnalysisTest, Construction ) {
	std::function<double( double )> dynamics;
	std::vector<hypro::Condition<double>> initialConditions;
	auto analyzer = analysis::Analyzer<hypro::Box<double>, decltype( dynamics )>( std::move( dynamics ), initialConditions );
	SUCCEED();
}