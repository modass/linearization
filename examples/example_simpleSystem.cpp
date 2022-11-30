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
 * Created by Stefan Schupp <stefan.schupp@tuwien.ac.at> on 18.11.22.
 */

#include "simple_example_matrices.h"
#include "simple_monomials.h"

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <functional>
#include <hypro/algorithms/reachability/Reach.h>
#include <hypro/datastructures/HybridAutomaton/HybridAutomaton.h>
#include <hypro/datastructures/HybridAutomaton/output/SpaceEx.h>
#include <hypro/datastructures/reachability/ReachTreev2.h>
#include <hypro/datastructures/reachability/Settings.h>
#include <hypro/flags.h>
#include <hypro/util/logging/Filewriter.h>
#include <hypro/util/plotting/Plotter.h>
#include <linearization/Conversion.h>
#include <linearization/mcCormick.h>
#include <vector>

int main( int argc, char** argv ) {
	// Command line interface
	double timeStepSize;
	double timeHorizon = 10;
	bool plot = false;
	std::size_t numberSystemVariables = 2;
	std::size_t numberSubdivisions = 5;

	CLI::App cli;
	cli.add_option( "-d, delta", timeStepSize, "Time step size used for reachability analysis" )->required()->check( CLI::PositiveNumber );
	cli.add_option( "-t, timeHorizon", timeHorizon, "Time horizon used for reachability analysis" )->check( CLI::PositiveNumber );
	cli.add_option( "-s, subdivisions", numberSubdivisions, "Number of grid points per dimension used to check the quality of the linearization when searching for a linearization point" )->check( CLI::PositiveNumber );
	cli.add_flag( "--plot", plot, "Set to create a plot of the set of reachable states projected to the first two dimensions" );
	try {
		cli.parse( argc, argv );
	} catch ( const CLI::ParseError& e ) {
		std::cout << cli.help() << std::endl;
		return cli.exit( e );
	}

	// using declarations
	using linearization::Domain;
	using linearization::Interval;
	using linearization::LinearizationResult;
	using linearization::MC;
	using linearization::Settings;

	// encode monomials as lambda functions
	auto monomialVector = linearization::order3Monomials<MC>();

	// settings for linearization
	auto intervalMonomials = linearization::order3Monomials<Interval>();
	std::vector<Interval> domains;
	std::vector<Interval> initialDomains;
	initialDomains.push_back( Interval{ 1, 1 } );
	initialDomains.push_back( Interval{ 2, 2 } );
	domains = initialDomains;
	// push dummy value
	for ( int i = numberSystemVariables; i < monomialVector.size(); ++i ) {
		initialDomains.push_back( Interval{ 0, 0 } );
	}
	std::transform( std::begin( intervalMonomials ), std::end( intervalMonomials ), std::back_inserter( domains ), [&initialDomains]( auto& f ) { return f( initialDomains ); } );
	Domain d{ domains };

	std::vector<std::size_t> subdivisions = std::vector<std::size_t>( numberSystemVariables + monomialVector.size() );
	std::fill_n( std::begin( subdivisions ), numberSystemVariables, numberSubdivisions );
	std::fill_n( std::next( std::begin( subdivisions ), numberSystemVariables ), monomialVector.size(), 1 );

	Settings s{ d, subdivisions };

	// combine initial constraints
	// add domain for x and y as well as a constraint system: constraints*x <= constants
	hypro::matrix_t<double> constraints = hypro::matrix_t<double>::Zero( 2 * numberSystemVariables, numberSystemVariables + monomialVector.size() );
	hypro::vector_t<double> constants = hypro::vector_t<double>::Zero( 2 * numberSystemVariables );
	for ( int i = 0; i < numberSystemVariables; ++i ) {
		constraints( 2 * i, i ) = 1;
		constants( 2 * i ) = d.intervals[i].u();
		constraints( 2 * i + 1, i ) = -1;
		constants( 2 * i + 1 ) = -d.intervals[i].l();
	}
	auto initialCondition = hypro::Condition<double>{ constraints, constants };

	// linearize
	for ( auto& monomial : monomialVector ) {
		LinearizationResult<double> linearizationResult;
		linearization::linearizeMonomial( monomial, s, linearizationResult );
		// add linearized initial sets for the other monomials
		initialCondition.addConstraints( linearizationResult.initialCondition );
	}

	spdlog::info( "Have computed over-approximation of the initial set: {}", initialCondition );

	// set up dynamic system as a hybrid automaton with one location only
	auto automaton = hypro::HybridAutomaton<double>();
	auto* loc = automaton.createLocation( "l0" );
	// encode dynamics: x' = Ax + b
	auto A = linearization::order3();
	loc->setFlow( hypro::linearFlow<double>( A.transpose() ) );

	automaton.setInitialStates( { { loc, initialCondition } } );

#ifdef HYPRO_ENABLE_SPACEEX_OUTPUT
	{
		hypro::LockedFileWriter xmlWriter{ "model.xml" };
		xmlWriter.clearFile();
		hypro::LockedFileWriter cfgWriter{ "config.cfg" };
		cfgWriter.clearFile();
		auto [xml, config] = hypro::toSpaceExFormat( automaton );
		xmlWriter << xml;
		cfgWriter << config;
	}
#endif

	spdlog::info( "Have created hybrid automaton: {}", automaton );

	// set up reachability analysis
	using Representation = hypro::SupportFunctionT<double, hypro::Converter<double>, hypro::NoBoxDetection>;
	auto roots = hypro::makeRoots<Representation, hypro::HybridAutomaton<double>>( automaton );
	// settings
	hypro::FixedAnalysisParameters fixedParameters{ 1, timeHorizon };
	hypro::AnalysisParameters parameters{ timeStepSize };
	hypro::reachability::Reach<Representation, hypro::HybridAutomaton<double>> reacher( automaton, fixedParameters, parameters, roots );

	// invoke analysis
	spdlog::info( "Start reachability analysis." );
	reacher.computeForwardReachability();

	// plot outputs
	if ( plot ) {
		spdlog::info( "Start plotting reachable sets." );
		auto& plt = hypro::Plotter<double>::getInstance();
		std::vector<std::size_t> plotDimensions{ 0, 1 };
		plt.setFilename( "simple_example" );
		plt.rSettings().overwriteFiles = false;
		for ( const auto& node : hypro::preorder( roots ) ) {
			for ( const auto& segment : node.getFlowpipe() ) {
				plt.addObject( segment.projectOn( plotDimensions ).vertices() );
			}
		}
		plt.plot2d( hypro::PLOTTYPE::png, true );
	}
	spdlog::info( "Computation finished, exiting." );
	return 0;
}
