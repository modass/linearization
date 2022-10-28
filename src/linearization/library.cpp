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

//supporting interval bounds are calculated using the default interval type mc::Interval here:
typedef mc::Interval I;
typedef mc::McCormick<I> MC;

//extremes of the intervals for the x and y variables.
double XL, XU, YL, YU;
//MSE for the concave and convex relaxation.
double sum1,sum2;

//definition of the function to be approximated.
template <class T>
T myfunc
( const T&x, const T&y )
{
  return 1.+x-sin(2.*x+3.*y)-cos(3.*x-5.*y);
}


//constructor.
MCpp::MCpp(double XLi, double XUi, double YLi, double YUi){
XL=XLi; XU=XUi; YL=YLi; YU=YUi;
NX=50;
NY=50;
std::cout <<XL<< '\n';
std::cout <<XU<< '\n';

}


// getters and setters
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

//MC++ library options
// provides tighter McCormick relaxations, but it is more time consuming.
MC::options.ENVEL_USE=true;
//maximum number of iterations for determination function points in convex/concave envelopes of univariate terms.
MC::options.ENVEL_MAXIT=100;
//tolerance for determination function points in convex/concave envelopes of univariate terms.
MC::options.ENVEL_TOL=1e-12;
//it indicates to use Tsoukalas & Mitsos's multivariate composition result for min/max, product, and division terms(more time consuming).
MC::options.MVCOMP_USE=true;

//output file
ofstream allsub( "MC-2D.out", ios_base::app);
allsub << scientific << setprecision(5) << right;

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------CONVEX RELAXATION  UNDERESTIMATING THE INITIAL NONLINEAR FUNCTION myfunc-----------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//we try 10 linearization points and then choose the one that minimizes the MSE between the affine approximation and McCormick relaxation at NX fixed points.
for (int i = 0; i < 10; ++i){
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    
/*set the point for which to compute the affine approximation (linearization point=(x_ref,y_ref)), sampling from a uniform distribution. The support for the considered uniform distribution is equal 
to the initial interval given for each coordinate.*/
    std::uniform_real_distribution<> disX(XL, XU);
    x_ref=disX(gen);
    std::uniform_real_distribution<> disY(YL, YU);
    y_ref=disY(gen);
    
    std::cout << x_ref << "x_ref_convex" << '\n';
    std::cout << y_ref << "y_ref_convex" << '\n';

    Xref=x_ref;
    Yref=y_ref;
    
    //Xrel and Yrel are variables of type McCormick, belonging to the interval [XL,XU] and [YL,YU] respectively, with current value is Xref and Yref.
    MC Xrel( I(XL,XU), Xref );
    MC Yrel( I(YL,YU), Yref );
    
    //Defining Xrel and Yrel as the subgradient components 0 and 1, respectively. 
    Xrel.sub(2,0);
    Yrel.sub(2,1);
    
    //McCormick relaxations of function myfunct on [XL,XU]\times[YL,YU] at (Xrel, Yrel):
    MC Zref = myfunc( Xrel, Yrel );
    cout << "Relaxation at reference point:\n" << Zref << endl;


   try{ 
   
    //compute the McCormick relaxations in NX points
    for( int iX=0; iX<NX; iX++ ){ 
      for( int iY=0; iY<NY; iY++ ){
        
        //set the coordinates of the points in which to compute the approximation, relying on NX and on the size of the initial intervals for each variable.
        double Xval = XL+iX*(XU-XL)/(NX-1.);
        double Yval = YL+iY*(YU-YL)/(NY-1.);
        double Zval = myfunc( Xval, Yval );
        
        //define Xrel and Yrel as variables of type McCormick, belonging to the interval [XL,XU] and [YL,YU] respectively, with current value is Xval and Yval. 
        MC Xrel( I(XL,XU), Xval );
        MC Yrel( I(YL,YU), Yval );
        
        //Defining Xrel and Yrel as the subgradient components 0 and 1, respectively. 
        Xrel.sub(2,0);
        Yrel.sub(2,1);
        //compute the McCormick convex and concave relaxations of myfunc at (Xrel,Yrel) along with subgradients of these relaxations.
        MC Zrel = myfunc( Xrel, Yrel );

       //The linearization point is chosen so that the distances between the MmcCormick relaxations and the affine approximation (computed in the sampled points) are minimized.
       //Here we update the partial sum for MSE computation, in particular we compute the difference between: 
       // - Zrel.cv(), that is the value of the McCormick convex relaxation at (Xval,Yval), and
       // - (Zref.cv()+Zref.cvsub(1)*(Yval-Yref)+Zref.cvsub(0)*(Xval-Xref)), that is the affine approximation generated using linearization point (Xref,Yref) and computed in (Xval,Yval)
        sum1 =sum1+pow(Zrel.cv()-(Zref.cv()+Zref.cvsub(1)*(Yval-Yref)+Zref.cvsub(0)*(Xval-Xref)),2);
      }}
    //takes the root to compute the MSE.
    sum1=sqrt(sum1);
    
    /*if the MSE obtained generating the affine approximation at the current linearization point leads to a lower MSE than the minimum obtained so far, then updates the final linearization point and the minimum MSE.*/
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
  
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------CONCAVE RELAXATION  OVERESTIMATING THE INITIAL NONLINEAR FUNCTION myfunc---------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//this code exactly mirrors the convex case  
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


//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------Here we save the results in the output file and compute the affine approximations corresponding to the concave and convex envelope for the chosen linearization points-------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  
 /*
 For each point at which we evaluate the function and the McCormick approximation, we also evaluate the affine approximation computed for the linearization point that minimizes the error made with respect to the mccormcick envelope
 */ 
 
  std::cout << "final:"<< Xref_c << '\n';
  std::cout << "final:"<< Yref_c << '\n';

  //compute the McCormick relaxation at the linearization point used for the affine underestimation.
  MC Xrel_c( I(XL,XU), Xref_c );
  MC Yrel_c( I(YL,YU), Yref_c );
  Xrel_c.sub(2,0);
  Yrel_c.sub(2,1);
  MC Zref_c = myfunc( Xrel_c, Yrel_c );
  cout << "Relaxation at reference point:\n" << Zref_c << endl;

  //compute the McCormick relaxation at the linearization point used for the affine overestimation.
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

        //print the coordinates of the points in which we are approximating the original function: Xval,Yval, the real value of the function at (Xval,Yval): Zval, the McCormick cancave and convex relaxation at  (Xval,Yval): Zrel.cc(),Zrel.cv(), and the values assumed by the affine over and under-approximation in (Xval,Yval).
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
  
