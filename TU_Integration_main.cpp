/*-------------------------------------------------------------------------------\
@ Numerical Methods by Yuha Park - Handong Global University

Author           : YuhaPark
Created          : 27-09-2024
Modified         : 03-10-2024
Language/ver     : C++ in MSVS2019

Description      : [Assignment] Integration_220002887_YuhaPark.cpp
-------------------------------------------------------------------------------*/

#include "stdio.h"
#include "stdlib.h"

#include "../../../include/myNP_22000287.h"

// Sample Code: Integration rectangular method 
double IntegrateRect(double x[], double y[], int m);

// You need to create myFunc() in this main source file
double myFunc(const double x);

void main()
{

	printf("\n**************************************************");
	printf("\n        PART 1. Integration from Datasets         ");
	printf("\n**************************************************\n");

	/************      Variables declaration & initialization      ************/
	double x[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };
	double y[] = { 0, 3, 8, 20, 33, 42, 40, 48, 60, 12, 8, 4, 3 };
	int M = sizeof(x) / sizeof(x[0]);


	/************      Solve  &	Show Output	   ************/
	double I_rect = IntegrateRect(x, y, M);
	printf("I_rect  = %f\n", I_rect);

	// Exercise 1. Trapezoid
	double I_trapz = 0;

	I_trapz = trapz(x, y, M);
	printf("I_trapz = %f\n\n", I_trapz);


	// Exercise 2. Simpson-------------------------------------------------------------

	/************      Variables declaration & initialization      ************/
	double x_[] = { -3, - 2.25, - 1.5, - 0.75,	0,	0.75,	1.5,	2.25,	3 };
	double y_[] = { 0,	2.1875,  3.75,	4.6875,  5,	4.6875,  3.75,	2.1875,	0 };
	int M_ = sizeof(x_) / sizeof(x_[0]);

	double I_simpson13 = 0;
	I_simpson13 = simpson13(x_, y_, M_);
	printf("I_simpson13  = %f\n\n", I_simpson13);



	printf("\n**************************************************");
	printf("\n        PART 2. Integration from a Function       ");
	printf("\n**************************************************\n");

	// Exercise 3. Integral
	double I_function = 0;

	I_function = integral(myFunc, -1, 1, 12);
	printf("I_function for Problem 3  = %f\n\n", I_function);

	I_function = integral(myFunc, -1, 1, 13); 
	printf("I_function for Problem 4  = %f\n\n", I_function);
	system("pause");
	
	
}



/*
y = sqrt(1-x^2)
*/
double myFunc(const double x) {

	double y = 0;

	y = sqrt(1 - x * x);

	return y;
}