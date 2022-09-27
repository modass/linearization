#include <iostream>
#include <iostream>
#include <random>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <string>
#include "boost/tuple/tuple.hpp"
#include "../../include/interval.hpp"
#include "../../include/mccormick.hpp"
#include "../linearization/include/library.h"

#define TEST_TRIG	
#define SAVE_RESULTS    
#undef USE_PROFIL	
#undef USE_FILIB


using namespace std;
using boost::tuple;
using namespace mc;

typedef mc::Interval I;
typedef mc::McCormick<I> MC;

double XL, XU, YL, YU;
double sum1,sum2;


template <class T>
T myfunc
( const T&x, const T&y )
{
  return 1.+x-sin(2.*x+3.*y)-cos(3.*x-5.*y);
}



MCpp::MCpp(double XLi, double XUi, double YLi, double YUi){
XL=XLi; XU=XUi; YL=YLi; YU=YUi;
NX=50;
NY=50;
std::cout <<XL<< '\n';
std::cout <<XU<< '\n';

}



void MCpp::setNX(int nx){
    NX = nx;
} 

void MCpp::setNY(int ny){ 
    NY = ny;
} 

int MCpp::getNX(){ 
    return NX;
} 

int MCpp::getNY(){ 
    return NX;
} 

void MCpp::setXBounds(int xl, int xu){ 
    XU = xu;
    XL = xl;
} 

void MCpp::setYBounds(int yl, int yu){ 
    YU = yu;
    YL = yl;
} 

int MCpp::getXBounds(){ 
    return XL,XU;
} 

int MCpp::getYBounds(){ 
    return YL,YU;
} 

int MCpp::do_relaxations() {
using std::cout;


double Xref, Yref, Xref_c, Yref_c, Xref_v, Yref_v;
double min=10000000000000000;
double x_ref = 0;
double y_ref = 0;

MC::options.ENVEL_USE=true;
MC::options.ENVEL_MAXIT=100;
MC::options.ENVEL_TOL=1e-12;
MC::options.MVCOMP_USE=true;


ofstream allsub( "MC-2D.out", ios_base::app);
allsub << scientific << setprecision(5) << right;


//--------------------------------CONVEX RELAXATION


for (int i = 0; i < 10; ++i){
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> disX(XL, XU);
    x_ref=disX(gen);
   
    std::uniform_real_distribution<> disY(YL, YU);
    y_ref=disY(gen);
    
    std::cout << x_ref << "x_ref_convex" << '\n';
    std::cout << y_ref << "y_ref_convex" << '\n';

    Xref=x_ref;
    Yref=y_ref;
    MC Xrel( I(XL,XU), Xref );
    MC Yrel( I(YL,YU), Yref );
    Xrel.sub(2,0);
    Yrel.sub(2,1);
    MC Zref = myfunc( Xrel, Yrel );
    cout << "Relaxation at reference point:\n" << Zref << endl;


   try{ 

    for( int iX=0; iX<NX; iX++ ){ 
      for( int iY=0; iY<NY; iY++ ){

        double Xval = XL+iX*(XU-XL)/(NX-1.);
        double Yval = YL+iY*(YU-YL)/(NY-1.);
        double Zval = myfunc( Xval, Yval );

        MC Xrel( I(XL,XU), Xval );
        MC Yrel( I(YL,YU), Yval );

        Xrel.sub(2,0);
        Yrel.sub(2,1);
        MC Zrel = myfunc( Xrel, Yrel );

        sum1 =sum1+pow(Zrel.cv()-(Zref.cv()+Zref.cvsub(1)*(Yval-Yref)+Zref.cvsub(0)*(Xval-Xref)),2);
      }}
    sum1=sqrt(sum1);
    if (sum1<min) {
    Xref_v=Xref;
    Yref_v=Yref;
    min=sum1;
}}
  
#ifndef USE_PROFIL
#ifndef USE_FILIB
  catch( I::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in natural interval extension:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }
#endif
#endif
  catch( MC::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in McCormick relaxation:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }
  sum1=0;
  }
  
//-------------------------------- CONCAVE RELAXATION
  
  sum2=0;
  min=10000000000000000;

for (int i = 0; i < 10; ++i){
  
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> disX(XL, XU);
    x_ref=disX(gen);
    std::uniform_real_distribution<> disY(YL, YU);
    y_ref=disY(gen);

    Xref=x_ref;
    Yref=y_ref;
 
    MC Xrel( I(XL,XU), Xref );
    MC Yrel( I(YL,YU), Yref );
    Xrel.sub(2,0);
    Yrel.sub(2,1);
    MC Zref = myfunc( Xrel, Yrel );
    cout << "Relaxation at reference point:\n" << Zref << endl;

   try{ 

    for( int iX=0; iX<NX; iX++ ){ 
      for( int iY=0; iY<NY; iY++ ){

        double Xval = XL+iX*(XU-XL)/(NX-1.);
        double Yval = YL+iY*(YU-YL)/(NY-1.);
        double Zval = myfunc( Xval, Yval );

        MC Xrel( I(XL,XU), Xval );
        MC Yrel( I(YL,YU), Yval );

        Xrel.sub(2,0);
        Yrel.sub(2,1);
        MC Zrel = myfunc( Xrel, Yrel );
        sum2=sum2+ pow((Zref.cc()+Zref.ccsub(1)*(Yval-Yref)+Zref.ccsub(0)*(Xval-Xref))-Zrel.cv(),2);      
      }

    }
    sum2=sqrt(sum2);
    if (sum2<min) {
    Xref_c=Xref;
    Yref_c=Yref;
    min=sum2;
}
  }
  
#ifndef USE_PROFIL
#ifndef USE_FILIB
  catch( I::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in natural interval extension:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }
#endif
#endif
  catch( MC::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in McCormick relaxation:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }
  sum2=0;
  }  

  //------------------------FINAL 
  std::cout << "final:"<< Xref_c << '\n';
  std::cout << "final:"<< Yref_c << '\n';

  MC Xrel_c( I(XL,XU), Xref_c );
  MC Yrel_c( I(YL,YU), Yref_c );
  Xrel_c.sub(2,0);
  Yrel_c.sub(2,1);
  MC Zref_c = myfunc( Xrel_c, Yrel_c );
  cout << "Relaxation at reference point:\n" << Zref_c << endl;
 
  MC Xrel_v( I(XL,XU), Xref_v );
  MC Yrel_v( I(YL,YU), Yref_v );
  Xrel_v.sub(2,0);
  Yrel_v.sub(2,1);
  MC Zref_v = myfunc( Xrel_v, Yrel_v );
  cout << "Relaxation at reference point:\n" << Zref_v << endl;
 

 try{ 

    for( int iX=0; iX<NX; iX++ ){ 
      for( int iY=0; iY<NY; iY++ ){

        double Xval = XL+iX*(XU-XL)/(NX-1.);
        double Yval = YL+iY*(YU-YL)/(NY-1.);
        double Zval = myfunc( Xval, Yval );

        MC Xrel( I(XL,XU), Xval );
        MC Yrel( I(YL,YU), Yval );

        Xrel.sub(2,0);
        Yrel.sub(2,1);
        MC Zrel = myfunc( Xrel, Yrel );

        allsub << setw(14) << Xval << setw(14) << Yval << setw(14) << Zval
               << setw(14) << Zrel.l() << setw(14) <<  Zrel.u()
               << setw(14) << Zrel.cv() << setw(14) << Zrel.cc()
               << setw(14) << Zref_v.cv()+Zref_v.cvsub(1)*(Yval-Yref_v)+Zref_v.cvsub(0)*(Xval-Xref_v)
               << setw(14) << Zref_c.cc()+Zref_c.ccsub(1)*(Yval-Yref_c)+Zref_c.ccsub(0)*(Xval-Xref_c)
               << endl;     
      }
      allsub << endl;
    }
  }
  
#ifndef USE_PROFIL
#ifndef USE_FILIB
  catch( I::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in natural interval extension:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }
#endif
#endif
  catch( MC::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in McCormick relaxation:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }


allsub.close();
return 0;
  }

  
