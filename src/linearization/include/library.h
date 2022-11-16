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

#ifndef LINEARIZATION_LIBRARY_H
#define LINEARIZATION_LIBRARY_H

class MCpp {
  public:
	MCpp( double XLi, double XUi, double YLi, double YUi );
	void setXBounds( int xl, int xu );
	void setYBounds(int yl, int yu);
int getXBounds(); 
int getYBounds();
void setNX(int nx);
void setNY(int ny);
int getNX();
int getNY();
int do_relaxations(); 

private:
double XL, XU, YL, YU;
int NX,NY;
double sum1,sum2;
  };
#endif//LINEARIZATION_LIBRARY_H

