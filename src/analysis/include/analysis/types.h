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

#ifndef ANALYSIS_TYPES_H
#define ANALYSIS_TYPES_H

#include <hypro/util/adaptions_carl/adaptions_includes.h>  // required for mpq_class and others

namespace analysis {

using Rational = mpq_class;
using Float = double;

}  // namespace analysis

#endif	// ANALYSIS_TYPES_H
