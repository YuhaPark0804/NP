/*-------------------------------------------------------------------------------\
@ C-Tutorial by Yuha-Park - Handong Global University
Author           : Yuha-Park 22000287
Created          : 08-30-2024
Modified         : 08-30-2024
Language/ver     : C++ in MSVS2022
Description      : TU_taylorSeries_exercise.cpp
-------------------------------------------------------------------------------*/


/* myNP_tutorial.cpp */

#include "myNP_tutorial.h"

void printVec(double* vec, int size)
{
	for (int i = 0; i < size; i++)
		printf("Vector[%d] = %.1f \n", i, vec[i]);
	printf("\n");
}


// factorial function
double factorial(int N)
{
	int y = 1;
	for (int k = 2; k <= N; k++)
		y = y * k;

	return y;
}

// power fuction
double power(double _x, int N)
{
	double y = 1;
	for (int i = 0; i < N; i++) {
		y *= _x; // y = x * y;
	}

	return y;
}


//  Taylor series approximation for sin(x) (input unit: [rad])
double sinTaylor(double _x)
{
	int N_max = 10;
	double S_N = 0;
	double S_k = 0;

	for (int k = 0; k < N_max; k++){
		S_k = power(-1, k) * power(_x, 2 * k + 1) / factorial(2 * k + 1);
		S_N += S_k;
	}

		return S_N;
}


// Taylor series approximation for sin(x) (input unit: [deg])
double sindTaylor(double _x)
{
	int N_max = 10;
	double S_N = 0;
	double S_k = 0;
	double x = (PI / 180) * _x;

	for (int k = 0; k < N_max; k++) {
		S_k = power(-1, k) * power(x, 2 * k + 1) / factorial(2 * k + 1);
		S_N += S_k;
	}

	return S_N;
}