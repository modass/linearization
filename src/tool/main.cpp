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

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <linearization/mcCormick.h>
#include <linearization/types.h>
#include <spdlog/spdlog.h>

int main( int argc, char** argv ) {
	using linearization::Domain;
	using linearization::Interval;
	using linearization::LinearizationResult;
	using linearization::MC;
	using linearization::Settings;

	// CLI
	CLI::App app{ "A tool for reachability analysis of non-linear dynamic systems via Carleman linearization." };
	// app.add_option("-f,--file", filename, "Help-string");
	CLI11_PARSE( app, argc, argv );

	// TODO perform Carleman linearization to obtain a set of monomials

	// convert to our data structures, which requires monomials to be a functor which takes a single vector as an input
	std::vector<std::function<MC( std::vector<MC> )>>
		  monomials;

	// perform linearization for each monomial
	for ( const auto& monomial : monomials ) {
		LinearizationResult<double> linearization;
		Settings linearization_settings;
		Domain d{ { Interval{ -2, 2 }, Interval{ 0, 4 } } };
		std::vector<std::size_t> subdivisions{ 5, 5 };
	}

	// construct linear system

	// run reachability analysis

	return 0;
}
