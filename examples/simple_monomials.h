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

#ifndef LINEARIZATION_SIMPLE_MONOMIALS_H
#define LINEARIZATION_SIMPLE_MONOMIALS_H

#include <functional>
#include <linearization/types.h>
#include <vector>

namespace linearization {

std::vector<std::function<MC( std::vector<MC> )>> order3Monomials() {
	// Matrix([[x], [y], [x**2*y], [x*y**2], [x**4*y], [x**3]])
	std::vector<std::function<MC( std::vector<MC> )>> result;
	result.push_back( []( const std::vector<MC>& in ) {MC res = pow( in[0], 2 ) * in[1] - in[2] + 0*in[3] + 0*in[4] + 0*in[5];
		return res; } );
	result.push_back( []( const std::vector<MC>& in ) {MC res = pow( in[1], 2 ) * in[0] - in[3] + 0*in[2] + 0*in[4] + 0*in[5];
		return res; } );
	result.push_back( []( const std::vector<MC>& in ) {MC res = pow( in[0], 4 ) * in[1] - in[4] + 0*in[2] + 0*in[3] + 0*in[5];
		return res; } );
	result.push_back( []( const std::vector<MC>& in ) {MC res = pow( in[0], 3 ) - in[5] + +0*in[1] + 0*in[2] + 0*in[3] + 0*in[4];
		return res; } );
	return result;
}

std::vector<std::function<MC( std::vector<MC> )>> order5Monomials() {
	// Matrix([[x], [y], [x**2*y], [x*y**2], [x**4*y], [x**3], [x**6*y], [x**3*y**2], [y**3], [x**5], [x**2*y**3], [x**5*y**2], [x**7], [x**8*y]])
	std::vector<std::function<MC( std::vector<MC> )>> result;
	result.push_back( []( const std::vector<MC>& in ) {MC res = pow( in[0], 2 ) * in[1] + 0*in[2] + 0*in[3] + 0*in[4] + 0*in[5] + 0*in[6] + 0*in[7] + 0*in[8] + 0*in[9] + 0*in[10] + 0*in[11] + 0*in[12] + 0*in[13];
		return res; } );
	result.push_back( []( const std::vector<MC>& in ) {MC res = in[0] * pow(in[1],2) + 0*in[2] + 0*in[3] + 0*in[4] + 0*in[5] + 0*in[6] + 0*in[7] + 0*in[8] + 0*in[9] + 0*in[10] + 0*in[11] + 0*in[12] + 0*in[13];
		return res; } );
	result.push_back( []( const std::vector<MC>& in ) {MC res = pow( in[0], 4 ) * in[1] + 0*in[2] + 0*in[3] + 0*in[4] + 0*in[5] + 0*in[6] + 0*in[7] + 0*in[8] + 0*in[9] + 0*in[10] + 0*in[11] + 0*in[12] + 0*in[13];
		return res; } );
	result.push_back( []( const std::vector<MC>& in ) {MC res = pow( in[0], 3 ) + 0*in[1] + 0*in[2] + 0*in[3] + 0*in[4] + 0*in[5] + 0*in[6] + 0*in[7] + 0*in[8] + 0*in[9] + 0*in[10] + 0*in[11] + 0*in[12] + 0*in[13];
		return res; } );
	result.push_back( []( const std::vector<MC>& in ) {MC res = pow( in[0], 6 ) * in[1] + 0*in[2] + 0*in[3] + 0*in[4] + 0*in[5] + 0*in[6] + 0*in[7] + 0*in[8] + 0*in[9] + 0*in[10] + 0*in[11] + 0*in[12] + 0*in[13];
		return res; } );
	result.push_back( []( const std::vector<MC>& in ) {MC res = pow( in[0], 3 ) * pow(in[1],2) + 0*in[2] + 0*in[3] + 0*in[4] + 0*in[5] + 0*in[6] + 0*in[7] + 0*in[8] + 0*in[9] + 0*in[10] + 0*in[11] + 0*in[12] + 0*in[13];
		return res; } );
	result.push_back( []( const std::vector<MC>& in ) {MC res = 0*in[0] + pow(in[1],3) + 0*in[2] + 0*in[3] + 0*in[4] + 0*in[5] + 0*in[6] + 0*in[7] + 0*in[8] + 0*in[9] + 0*in[10] + 0*in[11] + 0*in[12] + 0*in[13];
		return res; } );
	result.push_back( []( const std::vector<MC>& in ) {MC res = pow( in[0], 5 ) + 0*in[1] + 0*in[2] + 0*in[3] + 0*in[4] + 0*in[5] + 0*in[6] + 0*in[7] + 0*in[8] + 0*in[9] + 0*in[10] + 0*in[11] + 0*in[12] + 0*in[13];
		return res; } );
	result.push_back( []( const std::vector<MC>& in ) {MC res = pow( in[0], 2 ) * pow( in[1],3 ) + 0*in[2] + 0*in[3] + 0*in[4] + 0*in[5] + 0*in[6] + 0*in[7] + 0*in[8] + 0*in[9] + 0*in[10] + 0*in[11] + 0*in[12] + 0*in[13];
		return res; } );
	result.push_back( []( const std::vector<MC>& in ) {MC res = pow( in[0], 5 ) * pow( in[1],2 ) + 0*in[2] + 0*in[3] + 0*in[4] + 0*in[5] + 0*in[6] + 0*in[7] + 0*in[8] + 0*in[9] + 0*in[10] + 0*in[11] + 0*in[12] + 0*in[13];
		return res; } );
	result.push_back( []( const std::vector<MC>& in ) {MC res = pow( in[0], 7 ) + 0*in[1] + 0*in[2] + 0*in[3] + 0*in[4] + 0*in[5] + 0*in[6] + 0*in[7] + 0*in[8] + 0*in[9] + 0*in[10] + 0*in[11] + 0*in[12] + 0*in[13];
		return res; } );
	result.push_back( []( const std::vector<MC>& in ) {MC res = pow( in[0], 8 ) * in[1] + 0*in[2] + 0*in[3] + 0*in[4] + 0*in[5] + 0*in[6] + 0*in[7] + 0*in[8] + 0*in[9] + 0*in[10] + 0*in[11] + 0*in[12] + 0*in[13];
		return res; } );
	return result;
}

}  // namespace linearization

#endif	// LINEARIZATION_SIMPLE_MONOMIALS_H
