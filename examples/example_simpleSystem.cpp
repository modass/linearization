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

#include <functional>
#include <hypro/algorithms/reachability/Reach.h>
#include <hypro/datastructures/HybridAutomaton/HybridAutomaton.h>
#include <hypro/datastructures/reachability/ReachTreev2.h>
#include <hypro/datastructures/reachability/Settings.h>
#include <hypro/util/plotting/Plotter.h>
#include <linearization/Conversion.h>
#include <linearization/mcCormick.h>
#include <vector>

int main( int argc, char** argv ) {
	// using declarations
	using linearization::Domain;
	using linearization::Interval;
	using linearization::LinearizationResult;
	using linearization::MC;
	using linearization::Settings;
	// make readable names for variables
	const std::size_t x = 0;
	const std::size_t y = 1;
	const std::size_t z1 = 2;
	const std::size_t z2 = 3;
	const std::size_t z3 = 4;
	const std::size_t z4 = 5;

	// encode monomials as lambda functions
	/*
	Matrix([[x], [y], [x**2*y], [x*y**2], [x**4*y], [x**3]])
	*/
	std::function<MC( std::vector<MC> )> z1Monomial = []( const std::vector<MC>& in ) {MC res = pow( in[x], 2 ) * in[y] - in[z1] + 0*in[z2] + 0*in[z3] + 0*in[z4];
		return res; };
	std::function<MC( std::vector<MC> )> z2Monomial = []( const std::vector<MC>& in ) {MC res = pow( in[y], 2 ) * in[x] - in[z2] + 0*in[z1] + 0*in[z3] + 0*in[z4];
		return res; };
	std::function<MC( std::vector<MC> )> z3Monomial = []( const std::vector<MC>& in ) {MC res = pow( in[x], 4 ) * in[y] - in[z3] + 0*in[z1] + 0*in[z2] + 0*in[z4];
		return res; };
	std::function<MC( std::vector<MC> )> z4Monomial = []( const std::vector<MC>& in ) {MC res = pow( in[x], 3 ) - in[z4] + +0*in[y] + 0*in[z1] + 0*in[z2] + 0*in[z3];
		return res; };

	// settings for linearization
	auto xInterval = Interval{ 1, 1.5 };
	auto yInterval = Interval{ 2, 2.45 };
	auto z1Interval = pow( xInterval, 2 ) * yInterval;
	auto z2Interval = pow( yInterval, 2 ) * xInterval;
	auto z3Interval = pow( xInterval, 4 ) * yInterval;
	auto z4Interval = pow( xInterval, 3 );
	Domain d{ { xInterval, yInterval, z1Interval, z2Interval, z3Interval, z4Interval } };
	std::vector<std::size_t> subdivisions{ 5, 5, 1, 1, 1, 1 };
	Settings s{ d, subdivisions };

	// linearize
	LinearizationResult<double> z1LinearizationResult;
	linearization::linearizeMonomial( z1Monomial, s, z1LinearizationResult );
	LinearizationResult<double> z2LinearizationResult;
	linearization::linearizeMonomial( z2Monomial, s, z2LinearizationResult );
	LinearizationResult<double> z3LinearizationResult;
	linearization::linearizeMonomial( z3Monomial, s, z3LinearizationResult );
	LinearizationResult<double> z4LinearizationResult;
	linearization::linearizeMonomial( z4Monomial, s, z4LinearizationResult );

	// combine initial constraints
	auto initialCondition = hypro::Condition<double>{ z1LinearizationResult.initialCondition };
	initialCondition.addConstraints( z2LinearizationResult.initialCondition );
	initialCondition.addConstraints( z3LinearizationResult.initialCondition );
	initialCondition.addConstraints( z4LinearizationResult.initialCondition );
	// add domain for x and y as well as a constraint system: constraints*x <= constants
	{
		hypro::matrix_t<double> constraints = hypro::matrix_t<double>::Zero( 4, 6 );
		hypro::vector_t<double> constants = hypro::vector_t<double>::Zero( 4 );
		constraints( 0, x ) = 1;
		constants( 0 ) = xInterval.u();
		constraints( 1, x ) = -1;
		constants( 1 ) = -xInterval.l();
		constraints( 2, y ) = 1;
		constants( 2 ) = yInterval.u();
		constraints( 3, y ) = -1;
		constants( 3 ) = -yInterval.l();
		initialCondition.addConstraints( { constraints, constants } );
	}

	spdlog::info( "Have computed over-approximation of the initial set: {}", initialCondition );

	// set up dynamic system as a hybrid automaton with one location only
	auto automaton = hypro::HybridAutomaton<double>();
	auto* loc = automaton.createLocation( "l0" );
	// encode dynamics: x' = Ax + b
	/* A:
	[[ 0. -1.  0.  0.  0.  0.]
	[ 1.  1.  0.  0.  0.  0.]
	[ 0. -1.  1. -2.  0.  3.]
	[ 0.  0.  2.  2.  0.  0.]
	[ 0.  0. -1.  0.  1.  0.]
	[ 0.  0. -1.  0.  0.  0.]]
	 */
	hypro::matrix_t<double> A = hypro::matrix_t<double>::Zero( 7, 7 );
	A( 0, y ) = -1.0;
	A( 1, x ) = 1.0;
	A( 1, y ) = 1.0;
	A( 2, y ) = -1.0;
	A( 2, z1 ) = 1.0;
	A( 2, z2 ) = -2.0;
	A( 2, z4 ) = 3.0;
	A( 3, z1 ) = 2.0;
	A( 3, z2 ) = 2.0;
	A( 4, z1 ) = -1.0;
	A( 4, z3 ) = 1.0;
	A( 5, z1 ) = -1.0;
	loc->setFlow( hypro::linearFlow<double>( A ) );

	automaton.setInitialStates( { { loc, initialCondition } } );

	spdlog::info( "Have created hybrid automaton: {}", automaton );

	// set up reachability analysis
	using Representation = hypro::SupportFunctionT<double, hypro::Converter<double>, hypro::NoBoxDetection>;
	auto roots = hypro::makeRoots<Representation, double>( automaton );
	// settings
	hypro::FixedAnalysisParameters fixedParameters{ 1, 5 };
	hypro::AnalysisParameters parameters{ 0.01 };
	hypro::reachability::Reach<Representation> reacher( automaton, fixedParameters, parameters, roots );

	// invoke analysis
	reacher.computeForwardReachability();

	// plot outputs
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
	return 0;
}
