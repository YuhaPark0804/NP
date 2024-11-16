/* myNP_22000287.h */
#ifndef		_MY_NP_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NP_H
#define		PI		3.14159265358979323846264338327950288419716939937510582
#define M_NP  20
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


extern void printVec(double* vec, int row);


/*
calculates factorial
(_x)!
_x have to be int;
*/
extern double factorial(int _x);


/*
calculates power
_x^N;
N have to be int;
Therefore _x^0.5...(like sqrt) can not calculate.
*/
extern double power(double _x, int N);

/*returns sin(x) in unit[rad]
if x_0 is 0;
*/
extern double sinTaylor(double _x);

/*returns sin(x) in unit[deg]
if x_0 is 0;
*/
extern double sindTaylor(double _x);

//returns cos(x) in unit [rad]
double cosTaylor(double _x);

//returns cos(x) in unit [deg]
double cosdTaylor(double _x);

// Taylor series approximation for e^x
double exponential(double _x);


/* Bisection Method (solving Non Linear equations)
	_a      : initial value #1
	_b      : initial value #2
	_tol   : tolerance
*/
double bisection(double func(double _x), float _a0, float _b0, float _tol);

//double newtonRaphson(double _x0, double _tol);
double newtonRaphson(double func(double _x), double dfunc(double _x), double _x0, double _tol);

//secant_method
double secant_method(double func(double _x), double x0, double x1, double tol);

/*Differiantiation*/
// two-point difference approxiation: O(h) trucation error
void gradient1D(double x[], double y[], double dydx[], int m);

/*Differiantiation*/
// 1st order differentiation from the user-defined equation
void gradientFunc(double func(const double x), double _x[], double dydx[], int m);

/*Differiantiation
a function for 2nd order differentiation using discrete points

First Point: Four-point forward difference	(see APPENDIX)
Mid Points:	Three-point central difference
End Point:	Four-point backward difference
*/
void acceleration(double x[], double y[], double dy2dx2[], int m);

/*
편미분 1차
*/
void partialgradientFunc(double func(const double x, const double y), double x[], double y[], double dfdx[][M_NP], double dfdy[][M_NP], int m);
/*
편미분 2차
*/
void partial_2gradient_Func(double func(const double x, const double y), double x[], double y[], double d2fdx2[][M_NP], double d2fdy2[][M_NP], int m);

/*
Integration ----------------------------------------------------------------

Sample Code: Integration rectangular method

*/
double IntegrateRect(double x[], double y[], int m);

/*
integral by trapezoidal
*/
double trapz(double x[], double y[], int m);

/*
Integration Simpson 1/3 Method
Subinterval N = m - 1 have to be even

I = h/3 * [f(x0) + f(xN) + 4*sum(f(xi))_(i=1,3,5...N-1) + 2*sum(f(xk)_(k=2,4,6..N-2)]

*/
double simpson13(double x[], double y[], int m);

/*
Integration 
if (Subinterval N = m - 1 is even) Simpson 1/3 Method
if (Subinterval N = m - 1 is odd) Simpson 1/3 Method % 3/8 Method

3/8 Method: 

I = 3h/8 * [f(x0) + f(xN) + 3*sum(f(xi)+f(x_i+1)_(i=1,4,7...N-2) + 2*sum(f(xk)_(k=3,6,9..N-3)]
           
*/
double integral(double func(const double x), double a, double b, int n);


/*--------------------------ODE-------------------------*/
// 1st order ODE
void odeEU(double y[], double odeFunc(const double t, const double y), 
	const double t0, const double tf, const double h, const double y_init);
void odeRK2(double y[], double odeFunc(const double t, const double y), 
	const double t0, const double tf, const double h, const double y_init);
void odeRK3(double y[], double odeFunc(const double t, const double y), 
	const double t0, const double tf, const double h, const double y0);

// 2nd order ODE
void sys2RK2(double y[], double v[], void odeFuncSys(double dYdt[], const double t, const double Y[]), 
	const double t0, const double tf, const double h, const double y_init, const double v_init);
void sys2RK4(double y[], double v[], void odeFuncSys(double dYdt[], const double t, const double Y[]), 
	const double t0, const double tf, const double h, const double y_init, const double v_init);

/*------------------------End ODE-------------------------*/


#endif
