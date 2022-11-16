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

#include <iostream>
#include <random>
#include <stdlib.h>
#include "../linearization/include/library.h"


using namespace std;

int main() {
    /*
	double XL, XU, YL, YU;
    XL=-1;
    XU=1;  
    
    //Generation of a grid to subdivide the initial interval [XL,XU]\times[YL,YU] into subintervals to improve the accuracy of computed approximations.
    int Xspace=2;
    double Xstep=(XU-XL)/Xspace;
    cout << "-" << Xstep << "-" << endl;
    vector<double> arX;
    for (int i = 0; i < Xspace+1; i++) {
        arX.push_back(XL+i*Xstep);
        cout << arX[i] << " ";
    }
  
    YL=-2;
    YU=2;
    int Yspace=2;
    double Ystep=(YU-YL)/Yspace;
    cout << "-" << Ystep << "-" << endl;
    vector<double> arY;
    for (int i = 0; i < Yspace+1; i++) {
        arY.push_back(YL+i*Ystep);
        cout << arY[i] << " ";
    }
    
    //variables used to save the extremes of each subinterval.
    double xa,xb,ya,yb;
    
    //Compute the McCormick envelopes for each subinterval.
    for (int i = 1; i < arY.size(); i++) {
        ya=arY[i-1];
        yb=arY[i];
        for (int j = 1; j < arX.size(); j++) {
            xa=arX[j-1];
            xb=arX[j];
            MCpp solver(xa,xb,ya,yb);
            solver.do_relaxations();
    }
    }
     */

    return 0;
}
