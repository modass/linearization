#include <iostream>
#include <random>
#include <stdlib.h>
#include "../linearization/include/library.h"


using namespace std;

int main() {
    double XL, XU, YL, YU;	
    XL=-1;
    XU=1;   
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
    double xa,xb,ya,yb;
    
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

    return 0;
}





