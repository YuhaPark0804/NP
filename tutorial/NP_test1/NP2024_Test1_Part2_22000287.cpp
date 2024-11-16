/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [Yuha Park]
Created          : 10-07-2024
Modified         : 10-07-2024
Language/ver     : C++ in MSVS2022

Description      : MidTerm2024 - Part 2 (Q4)
---------------------------------------------------------------*/

// Include your library "myNP.h"
#include "stdio.h"
#include "stdlib.h"
#include <math.h>
#include "../../../include/myNP_22000287.h"


#define N 12

double func_Q4a(double x);
double func_Q4c(double x);

double integralQuad(double func(double _x));
double integralQuad2(double func(double _x), double a, double b);
double integralQuadComp(double func(double _x), double a, double b, int n);

int main(int argc, char* argv[])
{

	printf("\n**************************************************");
	printf("\n|                 Question 4-(a).                |");
	printf("\n**************************************************\n");

	double I1 = 0.0;

	double w0 = 5.0 / 9.0;
	double w1 = 8.0/ 9.0;
	double w2 = 5.0 / 9.0;

	double t0 = -1.0 * sqrt(1.0 / 3.0);
	double t1 = 0.0;
	double t2 = sqrt(1.0 / 3.0);
	

	I1 = integralQuad(func_Q4a);


	printf("\nanswer of Q4-(a) = %f\n\n", I1);




	printf("\n**************************************************");
	printf("\n|                 Question 4-(b).                |");
	printf("\n**************************************************\n");

	double a = 0.0;
	double b = 3.0;
	double I2 = 0.0;

	I2 = integralQuad2(func_Q4a, a, b);


	printf("\nanswer of Q4-(b) = %f\n\n", I2);



	printf("\n**************************************************");
	printf("\n|                 Question 4-(c).                |");
	printf("\n**************************************************\n");

	double A = 0.0;
	double B = 1.2;
	int    n = 12;
	double z_elong = 0.0;


	z_elong = integralQuadComp(func_Q4c, A, B, n);


	printf("\nanswer of Q4_c = %f\n\n", z_elong);


	system("pause");
	return 0;
}



/*==========================================================================*/
/*						Function Definitions								*/
/*==========================================================================*/


double func_Q4a(double x)
{

	double out = 0.0;
	
	out = exp(-1*pow(x, 2));

	return out;
}

double func_Q4c(double x)
{
	double out = 0.0;

	double ep = 40 * pow(10, -3) * pow(x, 0.5);
	double L = 1.2; //[m]

	out = (1 + ep);

	return out;
}

double integralQuad(double func(double _x))
{
	double weight[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
	double xn[3] = { -sqrt(1.0 / 3.0), 0.0, sqrt(1.0 / 3.0) };

	double output = 0.0;
	double output1 = 0.0;
	double output2 = 0.0;
	double output3 = 0.0;

	for (int i = 0; i < 3; i++) {
		output += weight[i] * func(xn[i]);
	}

	return output;
}

double integralQuad2(double func(double _x), double a, double b)
{
	double weight[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
	double xn[3] = { -sqrt(1.0 / 3.0), 0.0, sqrt(1.0 / 3.0) };

	double output = 0.0;


	for (int i = a; i < b; i++) {
		output += weight[i] * func(0.5 * ((b - a) * xn[i] + a + b)) * 0.5 * (b - a);
	}
	
	return output;
}

double integralQuadComp(double func(double _x), double a, double b, int n)
{

	double weight[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
	double xn[3] = { -sqrt(1.0 / 3.0), 0.0, sqrt(1.0 / 3.0) };

	double output[N] = { 0.0 };
	double out = 0.0;

	for (int k = 0; k < n; k++) {
		for (int i = a; i < b; i++) {
			output[k] += weight[i] * func(0.5 * ((b - a) * xn[i] + a + b)) * 0.5 * (b - a);
		}
	}

	for (int i = 0; i < n; i++) {
		out += output[i];
	}


	return out;

}
