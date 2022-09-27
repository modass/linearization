#ifndef LINEARIZATION_LIBRARY_H
#define LINEARIZATION_LIBRARY_H

class MCpp
{

public:
MCpp(double XLi, double XUi, double YLi, double YUi);
void setXBounds(int xl, int xu);
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

