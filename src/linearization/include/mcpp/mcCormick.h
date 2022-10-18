/*
 * Copyright (c) 2022 Stefan Schupp.
 * This file is part of the linearization project.
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
 *
 */

//
// Created by Stefan Schupp on 17.10.22.
//

#ifndef LINEARIZATION_MCCORMICK_H
#define LINEARIZATION_MCCORMICK_H

#include <functional>
#include <hypro/representations/GeometricObjectBase.h>
#include <utility/Exceptions.h>

namespace linearization
{

/*
 * Class that wraps mcpp
 */
class MCCormick
{
public:
    explicit MCCormick(std::function<double()>&& dynamics)
        : mDynamics(std::move(dynamics))
    {
    }

    /**
     * Computes a polytope over-approximating a function within a given interval
     * @return The polytope
     */
    hypro::HPolytope<double> linearize() const;

private:
    std::function<double()> mDynamics;
};

} // namespace linearization

#endif // LINEARIZATION_MCCORMICK_H
