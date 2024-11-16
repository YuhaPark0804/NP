/* myNP_tutorial.h */
#ifndef		_MY_NP_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NP_H
#define		PI		3.14159265358979323846264338327950288419716939937510582

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


extern void printVec(double* vec, int row);


//calculates factorial
extern double factorial(int _x);

// calculates power
extern double power(double _x, int N);

//returns sin(x) in unit [rad]
extern double sinTaylor(double _x);

//returns sin(x) in unit [deg]
extern double sindTaylor(double _x);

#endif
