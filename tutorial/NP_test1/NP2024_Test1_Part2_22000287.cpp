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
#include "../../include/myNP_22000287.h"


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

	// [TO-DO] Add your code here
	// [TO-DO] Add your code here
	// [TO-DO] Add your code here

	printf("\nanswer of Q4-(a) = %f\n\n", I1);



	
	printf("\n**************************************************");
	printf("\n|                 Question 4-(b).                |");
	printf("\n**************************************************\n");
	
	double a    = 0.0;
	double b    = 3.0;
	double I2 = 0.0;

	// [TO-DO] Add your code here
	// [TO-DO] Add your code here
	// [TO-DO] Add your code here

	printf("\nanswer of Q4-(b) = %f\n\n", I2);
	
	


	
	printf("\n**************************************************");
	printf("\n|                 Question 4-(c).                |");
	printf("\n**************************************************\n");
	
	double a    = 0.0;
	double b    = 1.2;
	int    n    = 12;
	double z_elong = 0.0;

	
	// [TO-DO] Add your code here
	// [TO-DO] Add your code here
	// [TO-DO] Add your code here
	
		
	printf("\nanswer of Q4_c = %f\n\n", z_elong);


	system("pause");
	return 0;
}



/*==========================================================================*/
/*						Function Definitions								*/
/*==========================================================================*/


double func_Q4a(double x)
{
	
	int N_MAX = 10;
	double out = 0.0;
	double x = _x;

	for (int k = 0; k < N_MAX; k++) {
		out += pow(x, k) / factorial(k);
	}
	return S_N;

	return out;
}

double func_Q4c(double x)
{
	double out = 0.0;

	// [TO-DO] Add your code here
	// [TO-DO] Add your code here

	return out;
}

double integralQuad(double func(double _x))
{
	double weight[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
	double xn[3]     = {-sqrt(1.0 / 3.0), 0.0, sqrt(1.0 / 3.0)};

	double output = 0.0;

	// [TO-DO] Add your code here
	// [TO-DO] Add your code here

	return output;
}

double integralQuad2(double func(double _x), double a, double b)
{
	//double weight[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
	//double xn[3] = { -sqrt(1.0 / 3.0), 0.0, sqrt(1.0 / 3.0) };

	double output = 0.0;

	// [TO-DO] Add your code here
	// [TO-DO] Add your code here

	return output;
}

double integralQuadComp(double func(double _x), double a, double b, int n)
{
	double output = 0.0;

	// [TO-DO] Add your code here
	// [TO-DO] Add your code here

	return output;

}
