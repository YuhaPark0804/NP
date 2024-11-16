/*-------------------------------------------------------------------------------\
@ C-Tutorial by Yuha-Park - Handong Global University
Author           : Yuha-Park 22000287
Created          : 30-08-2024
Modified         : 14-11-2024
Language/ver     : C++ in MSVS2022
-------------------------------------------------------------------------------*/


/* myNP_220002878.cpp */

#include "myNP_22000287.h"

void printVec(double* vec, int size) //--------------------------------------------------
{
	for (int i = 0; i < size; i++)
		printf("Vector[%d] = %.1f \n", i, vec[i]);
	printf("\n");
}


// factorial function----------------------------------------------------------------------
double factorial(int N)
{
	int y = 1;
	for (int k = 2; k <= N; k++)
		y = y * k;

	return y;
}

// power fuction------------------------------------------------------------------------------
double power(double _x, int N)
{
	double y = 1;
	for (int i = 0; i < N; i++) {
		y *= _x; // y = x * y;
	}

	return y;
}


//  Taylor series approximation for sin(x) (input unit: [rad])----------------------------------
double sinTaylor(double _x)
{
	int N_max = 10;
	double S_N = 0;
	double S_k = 0;

	for (int k = 0; k < N_max; k++) {
		S_k = power(-1, k) * power(_x, 2 * k + 1) / factorial(2 * k + 1);
		S_N += S_k;
	}

	return S_N;
}


// Taylor series approximation for sin(x) (input unit: [deg])-------------------------------
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

//  Taylor series approximation for cos(x) (input unit: [rad])----------------------------------
double cosTaylor(double _x)
{
	int N_max = 10;
	double S_N = 0;
	double S_k = 0;

	for (int k = 0; k < N_max; k++) {
		S_k = power(-1, k) * power(_x, 2 * k) / factorial(2 * k);
		S_N += S_k;
	}

	return S_N;
}

// Taylor series approximation for cos(x) (input unit: [deg])-------------------------------------
double cosdTaylor(double _x)
{
	int N_max = 10;
	double S_N = 0;
	double S_k = 0;
	double x = (PI / 180) * _x;

	for (int k = 0; k < N_max; k++) {
		S_k = power(-1, k) * power(_x, 2 * k) / factorial(2 * k);
		S_N += S_k;
	}

	return S_N;
}

// Taylor series approximation for e^x--------------------------------------------------------------
double exponential(double _x) {
	int N_MAX = 10;
	double S_N = 0;
	double x = _x;

	for (int k = 0; k < N_MAX; k++) {
		S_N += pow(x, k) / factorial(k);
	}
	return S_N;
}






/* Bisection Method-------------------------------------------------------------------------------
	_a      : initial value #1
	_b      : initial value #2
	_tol   : tolerance

	ex)
	printf("Bisection Method:\n");
	sol_bm = bisection(func, a0, b0, tol);

	printf("Final Solution: %f \t", sol_bm);;
*/
double bisection(double func(double _x), float _a0, float _b0, float _tol)
{
	// Define Bisection function, assuming (func(a) * func(b) <0 )
	// Initialization
	int k = 0;
	int Nmax = 1000;
	float a = _a0;
	float b = _b0;
	float xn = 0;
	float ep_old = 1000;
	float ep = 1000;
	int count = 0;

	// Bisection 
	while (k<Nmax && ep>_tol) {
		// Update xn as midpoint
		xn = 0.5 * (a + b);

		// Update range a, b
		if (func(xn) * func(a) < 0) {
			b = xn;
		}
		else if (func(xn) * func(b) < 0){
			a = xn;
		}
		

		// Check tolerance
		ep = fabs(func(xn));
		if (ep == ep_old){
			if (count == 1) break;
			count++;
		}
		ep_old = ep;

		k++;

		printf("k:%d \t", k);
		printf("Xn(k): %f \t", xn);
		printf("Tol: %.10f\n", ep);
	}

	return xn;
}



/*
Newton-Raphson Method without Passing a Function ----------------------------------------------------- 

ex)
	printf("Newton-Raphson Method Result:\n");
	sol_nr = newtonRaphson(func, dfunc, x0, tol);

	printf("Final Solution: %f \t", sol_nr);
*/
double newtonRaphson(double func(double _x), double dfunc(double _x), double _x0, double _tol)
{
	//float xn = _x0;
	float ep = 1000;
	int Nmax = 1000;
	int k = 0;
	double h = 0;

	float x_pre = 0;
	float  x_old = _x0;
	int count = 0;
	float ep_old = 100;

	while (k<Nmax && ep>_tol) {
		if (dfunc(x_old) == 0)
		{
			printf("[ERROR] dF == 0 !!\n");
			break;
		}
		else
		{
			                                       // get h=f/df @ x(k)
			h = - func(x_old) / dfunc(x_old);

			                                       // update x(k+1)=x(k)+h(k)
			//xn += h;
			x_pre = x_old + h;

												     // check tolerance
			ep = fabs(func(x_pre));		             //ep check method 1
			//ep = fabs((x_pre - x_old) / x_old);    //ep check method 2
			
			if (ep == ep_old) {
				if (count == 1) break;
				count++;
			}
			x_old = x_pre;
			ep_old = ep;

			k++;

			printf("k:%d \t", k);
			printf("X(k): %f \t", x_pre);
			printf("Tol: %.10f\n", ep);

		}
	}

	return x_pre;
}

/*Secant method----------------------------------------------------------------------------

ex)
	printf("Secant Method Result:\n");
	sol_sm = secant_method(func, Ans0, Ans1, tol);

	printf("Final Solution: %f \t", sol_sm);
*/

double secant_method(double func(double _x), double x0, double x1, double tol) {
	int NMAX = 1000;
	double ep = 1000;
	int k = 0;

	double x_pre = 0;
	double x_now = x1;
	double x_old = x0;

	double ep_old = 0;
	int count = 0;

	while (k < NMAX && ep > tol) {
		
		if (func(x_now) - func(x_old) == 0) {
			printf("[ERROR] func(x_now) - func(x_old) == 0 !!\n");
			break;
		}
		
		x_pre = x_now - func(x_now) * (x_now - x_old) / (func(x_now) - func(x_old));

		ep = fabs(func(x_pre));

		if (ep == ep_old) {
			if (count == 1) break;
			count++;
		}
		ep_old = ep;

		x_old = x_now;
		x_now = x_pre;
		
		k++;

		printf("k:%d \t", k);
		printf("X(k): %f \t", x_pre);
		printf("Tol: %.10f\n", ep);

	}

	return x_pre;
}


/*two-point difference approxiation: O(h) trucation error------------------------------------------------

ex)
    int m = 12;
	double t[12] = { 0 };
	for (int i = 0; i < m; i++) t[i] = 0.5 * (i - 2);

	double x[] = { -3.632, -0.3935, 1, 0.6487, -1.282, -4.518, -8.611, -12.82, -15.91, -15.88, -9.402, 9.017 };

	double  dxdt[12] = { 0 };

	gradient1D(t, x, dxdt, m);
	printVec(dxdt, m);

*/

void gradient1D(double x[], double y[], double dydx[], int m) {

	double h = x[1] - x[0];

	dydx[0] = (-3 * y[0] + 4 * y[1] - y[2]) / (2 * h);                          // for the first point x[0]; Three point forward difference

	for (int i = 1; i < m - 1; i++){                                          
		dydx[i] = (y[i + 1] - y[i - 1]) / ((2 * h));                      // for x[1]~x[m-2]; Two point central difference
	}

	dydx[m - 1] = (y[m - 3] - 4 * y[m - 2] + 3 * y[m - 1]) / ((2 * h));   // for the last point x[m-1]; Three point backward difference
}



/*
1st order differentiation from the user-defined equation---------------------------------------
	#define M 10

	double h = 0.01;
	double x_Q2[M] = { 2.0, };

	for (int i = 1; i < M; i++) { 
		x_Q2[i] = x_Q2[i - 1] + h;
	}

	double dydx_Q2 = 0;
	double dydx[M] = { 0.0 };

	gradientFunc(func_Q2, x_Q2, dydx, M);
	dydx_Q2 = dydx[0];

	printf("\nanswer of Q2 dydx= %f\n\n", dydx_Q2);
	
*/

/*
double xin = 2;
	m = 21;
	double dydx[21] = { 0 };  // m = 21 points
	double dy2dx2[21] = { 0 };  // m = 21 points

	// User defined function F(x)
	double y = myFunc(xin);
	printf("\n y=myFun(xin) = %f \n\n", y);

	double _t[21] = { 0 };
	for (int i = 0; i < 21; i++) _t[i] = 0.2 * i;  // 0,0.2,0.4 ... 4

	gradientFunc(myFunc, _t, dydx, m);
	printVec(dydx, m);
*/


void gradientFunc(double func(const double x), double x[], double dydx[], int m) {

	int M = m;
	
	dydx[0] = (-3 * func(x[0]) + 4 * func(x[1]) - func(x[2])) / (x[2] - x[0]);                          // for the first point x[0]

	for (int i = 1; i < m - 1; i++) {
		dydx[i] = (func(x[i + 1]) - func(x[i - 1])) / (x[i + 1] - x[i - 1]);                            // for x[1]~x[m-2]
	}

	dydx[m - 1] = (func(x[m - 3]) - 4 * func(x[m - 2]) + 3 * func(x[m - 1])) / (x[m - 1] - x[m - 3]);   // for the last point x[m-1]

	
}


/*Differiantiation------------------------------------------------------------------
a function for 2nd order differentiation using discrete points

First Point: Four-point forward difference	(see APPENDIX)
Mid Points:	Three-point central difference
End Point:	Four-point backward difference
*/

/* ex)

  double _y[21] = { 0 };
	for (int i = 0; i < 21; i++) {
		_y[i] = myFunc(_t[i]);
	}

	acceleration(_t, _y, dy2dx2, m);
	printVec(dy2dx2, m);
*/
void acceleration(double x[], double y[], double dy2dx2[], int m) {

	dy2dx2[0] = (2 * y[0] - 5 * y[1] + 4 * y[2] - y[3]) / pow((x[1] - x[0]), 2);                                // for the first point x[0]

	for (int i = 1; i < m - 1; i++) {
		dy2dx2[i] = (y[i - 1] - 2 * y[i] + y[i + 1]) / pow((x[i + 1] - x[i]), 2);                                     // for x[1]~x[m-2]
	}

	dy2dx2[m - 1] = (-y[m - 4] + 4 * y[m - 3] - 5 * y[m - 2] + 2 * y[m - 1]) / pow((x[m - 1] - x[m - 2]), 2);   // for the last point x[m-1]
}

/*
 Æí¹ÌºÐ

	
	double dydx_t[M_NP][M_NP] = { 0 };
	double dydt_t[M_NP][M_NP] = { 0 };
	double dy2dx2[M_NP][M_NP] = { 0 };
	double dy2dt2[M_NP][M_NP] = { 0 };
	double t[M_NP] = { 0 }; //1
	double x[M_NP] = { 0 }; //0.2


	double a_t = 0.0, b_t = 2.0;
	double a_x = 0.0, b_x = 1.0;
	double h_t = (b_t - a_t) / (M_NP); //0.1
	double h_x = (b_x - a_x) / (M_NP); //0.05

	for (int i = 0; i < M_NP; i++) {  // 0, 0.1,...,1.9, 2.0;   t[10] = 1.0
		t[i] = a_t + i * h_t;
		printf("t[%d]: %f\n", i, t[i]);

	}
	for (int i = 0; i < M_NP; i ++) {  // 0, 0.05, ...0.95, 1.0;  x[4] = 0.2
		x[i] = a_x + i * h_x;
		printf("x[%d]: %f\n", i, x[i]);
	}


	partialgradientFunc(myfunc4, x, t, dydx_t, dydt_t, M_NP);
	partial_2gradient_Func(myfunc4, x, t, dy2dx2, dy2dt2, M_NP);


	printf("\nanswer of Q3 a) (t=1, x=0.2) = %f\t%f\n", dydt_t[4][10], dydx_t[4][10]);
	printf("\nanswer of Q3 b) (t=1, x=0.2) = %f\t%f\n", dy2dt2[4][10], dy2dx2[4][10]);
	
*/
void partialgradientFunc(double func(const double x, const double y), double x[], double y[], double dfdx[][M_NP], double dfdy[][M_NP], int m) {
	
	double hx = x[1] - x[0]; //0.05
	double hy = y[1] - y[0]; //0.1

	// dfdx[]

	for (int j = 0; j < m - 1; j++) {
		dfdx[0][j] = (func(x[1], y[j]) - func(x[0], y[j])) / (hx);  //forward
	}

	for (int i = 1; i < m - 1; i++) {
		for (int j = 0; j < m - 1; j++) {
			dfdx[i][j] = (func(x[i + 1], y[j]) - func(x[i - 1], y[j])) / (2 * hx);  //central
		}
	}
	

	for (int j = 0; j < m - 1; j++) {
		dfdx[m - 1][j] = (func(x[m - 1], y[j]) - func(x[m - 2], y[j])) / (hx);  //backward
	}

	// dfdy[]

	for (int i = 0; i < m - 1; i++) {
		dfdy[i][0] = (func(x[i], y[1]) - func(x[i], y[0])) / (hy);
	}

	for (int i = 0; i < m - 1; i++) {
		for (int j = 1; j < m - 1; j++) {
			dfdy[i][j] = (func(x[i], y[j + 1]) - func(x[i], y[j - 1])) / (2 * hy);
		}

	}
	for (int i = 0; i < m - 1; i++) {
		dfdy[i][m - 1] = (func(x[i], y[m - 1]) - func(x[i], y[m - 2])) / (hy);
	}

}


/*
ex) central

for (int i = 1; i < m - 1; i++) {
		for (int j = 0; j < m - 1; j++) {
			d2fdxdy[i][j] = (func(x[i + 1], t[j + 1]) - func(x[i - 1], t[j + 1]) - (func(x[i + 1], t[j - 1]) - func(x[i - 1], t[i - 1]))) / (2* hx * 2 * hy);
		}
	}

*/
void partial_2gradient_Func(double func(const double x, const double y), double x[], double t[], double d2fdx2[][M_NP], double d2fdy2[][M_NP], int m) {
	double hx = x[1] - x[0];
	double hy = t[1] - t[0];

	//dy2dx2

	for (int i = 1; i < m - 1; i++) {
		for (int j = 0; j < m - 1; j++) {
			d2fdx2[i][j] = (func(x[i - 1], t[j]) - 2 * func(x[i], t[j]) + func(x[i + 1], t[j])) / (hx * hx);
		}
	}

	//dy2dt2

	for (int i = 0; i < m - 1; i++) {
		for (int j = 1; j < m - 1; j++) {
			d2fdy2[i][j] = (func(x[i], t[j - 1]) - 2 * func(x[i], t[j]) + func(x[i], t[j + 1])) / (hy * hy);
		}
		
	}

}








/*
Integration ----------------------------------------------------------------

Sample Code: Integration rectangular method 

*/
/*
	double x[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };
	double y[] = { 0, 3, 8, 20, 33, 42, 40, 48, 60, 12, 8, 4, 3 };
	int M = sizeof(x) / sizeof(x[0]);

	double I_rect = IntegrateRect(x, y, M);
	printf("I_rect  = %f\n", I_rect);
*/
 
double IntegrateRect(double x[], double y[], int m) {
	int N = m - 1;
	double I = 0;
	for (int i = 0; i < N; i++)
		I += y[i] * (x[i + 1] - x[i]);
	return I;
}

/*
Integration ----------------------------------------------------------------

Trapezoidal Method
assume x[i+1] - x[i] = h (const);
}*/

/*
	double x[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };
	double y[] = { 0, 3, 8, 20, 33, 42, 40, 48, 60, 12, 8, 4, 3 };
	int M = sizeof(x) / sizeof(x[0]);

	double I_trapz = 0;

	I_trapz = trapz(x, y, M);
	printf("I_trapz = %f\n\n", I_trapz);
*/
double trapz(double x[], double y[], int m) {
	double h = x[1] - x[0];
	double F = 0;
	double I = 0;
	int N = m - 1;        //N: intervals, m: data points

	for (int idx = 1; idx < N ; idx++) {
		I += y[idx];
	}
	F = 0.5 * h * (y[0] + y[N]) + h * I;  // h/2*(y[0] + y[N]) + h*sum(y[1~N-1])

	return F;

}


/*
Integration Simpson 1/3 Method
Subinterval N = m - 1 have to be even!!!!!!!!

I = h/3 * [f(x0) + f(xN) + 4*sum(f(xi))_(i=1,3,5...N-1) + 2*sum(f(xk)_(k=2,4,6..N-2)]
*/

/*
	double x_[] = { -3, - 2.25, - 1.5, - 0.75,	0,	0.75,	1.5,	2.25,	3 };
	double y_[] = { 0,	2.1875,  3.75,	4.6875,  5,	4.6875,  3.75,	2.1875,	0 };
	int M_ = sizeof(x_) / sizeof(x_[0]);  //9

	double I_simpson13 = 0;
	I_simpson13 = simpson13(x_, y_, M_);
	printf("I_simpson13  = %f\n\n", I_simpson13);
*/
double simpson13(double x[], double y[], int m) {
	double h = x[1] - x[0];
	double I1 = 0;
	double I2 = 0;
	double I = 0;
	int N = m - 1;

	
	for (int i = 1; i < N ; i += 2) {  //sum(f(xi))_(i=1,3,5...N-1)
		I1 += y[i];
	}

	for (int j = 2; j < N - 1; j += 2) {  // sum(f(xk)_(k=2,4,6..N-2)
		I2 += y[j];
	}

	I = y[0] + y[N] + (4.0 * I1) + (2.0 * I2);
	I *= (h / 3.0);

	return I;
}
	
/*
Integration-----------------------------------------------------------------------------
if (Subinterval N = m - 1 is even) Simpson 1/3 Method
if (Subinterval N = m - 1 is odd) Simpson 1/3 Method % 3/8 Method

h = (b - a) / n (h is const) (smae intervals)

1/3 Method:
I = h/3 * [f(x0) + f(xN) + 4*sum(f(xi))_(i=1,3,5...N-1) + 2*sum(f(xk)_(k=2,4,6..N-2)]

3/8 Method:
I = 3h/8 * [f(x0) + f(xN) + 3*sum(f(xi)+f(x_i+1)_(i=1,4,7...N-2) + 2*sum(f(xk)_(k=3,6,9..N-3)]

*/

/*
	double I_function = 0;

	I_function = integral(myFunc, -1, 1, 12);
	printf("I_function for Problem 3  = %f\n\n", I_function);

	I_function = integral(myFunc, -1, 1, 13);
	printf("I_function for Problem 4  = %f\n\n", I_function);
*/
double integral(double func(const double x), double a, double b, int n) {

	double h = (b - a) / n;
	double I1 = 0;
	double I2 = 0;
	double I = 0;
	int idx = 0;
	int N = n - 3; // if(n=13) 10 = 13 - 3

	double* x, * y;
	x = (double*)malloc(sizeof(double) * (n + 1));  //if n=12, array: 0~12. so, size = n + 1
	y = (double*)malloc(sizeof(double) * (n + 1));

	for (double i = a; i <= b; i += h) {
		x[idx] = i;             
		y[idx] = func(i);        
		idx++;
	}


	if (n % 2 == 0) {  // if Subinterval(n) is even; use Simpson 1/3 Method only
		for (int i = 1; i < n; i += 2) {      //sum(f(xi))_(i=1,3,5...N-1)
			I1 += y[i];
		}

		for (int j = 2; j < n - 1; j += 2) {  //sum(f(xk)_(k=2,4,6..N-2)
			I2 += y[j];
		}

		I = y[0] + y[n] + (4.0 * I1) + (2.0 * I2);
		I *= (h / 3.0);

		free(x);
		free(y);

		return I;
	}
	else {         // if Subinterval(n) is odd; use Simpson 1/3 & 3/8 Method 
		for (int i = 1; i < N; i += 2) {  //sum(f(xi))_(i=1,3,5...N-1)
			I1 += y[i];
		}

		for (int j = 2; j < N - 1; j += 2) {      //sum(f(xk)_(k=2,4,6..N-2)
			I2 += y[j];
		}

		I += (h / 3.0) * (y[0] + y[N] + (4.0 * I1) + (2.0 * I2));
		I += (3.0 / 8.0) * h * (y[N] + 3.0 * y[N + 1] + 3.0 * y[N + 2] + y[N + 3]);

		free(x);
		free(y);

		return I;
	}
}

/*------------------ODE--------------------*/

void odeEU(double y[], double odeFunc(const double t, const double y), const double t0, const double tf, const double h, const double y_init) {

	/*Input:
	*	y[]:			A vector with the y coordinate of the solution points. 
	*	odeFunc(t, y):	Name of a function file that calculates dy/dt.  
	*	t0:				The first value of t.
	*	tf:				The last value of t.
	*	h:				Step size.
	*	y_init:			Initial Condition of y(1)
	*/
	
	// Variable Initialization
	int N = (tf - t0) / h + 1;

	//Initial Condition
	double ti = t0;
	y[0] = y_init;

	// Euler Explicit ODE Method
	for (int i = 0; i < N - 1; i++) {
		ti += h;
		y[i + 1] = y[i] + odeFunc(ti, y[i]) * h;
	}
}

void odeRK2(double y[], double odeFunc(const double t, const double y), const double t0, const double tf, const double h, const double y_init) {

	/*Input:
	*	y[]:			A vector with the y coordinate of the solution points.
	*	odeFunc(t, y):	Name of a function file that calculates dy/dt.
	*	t0:				The first value of t.
	*	tf:				The last value of t.
	*	h:				Step size.
	*	y_init:			Initial Condition of y(1)
	*/

	// Variable Initialization
	int N = (tf - t0) / h + 1;
	double y2 = 0;
	double t2 = 0;

	// Initial Condition
	double ti = t0;
	y[0] = y_init;

	// RK Design Parameters
	double C1 = 0.5;
	double C2 = 1 - C1;
	double alpha = 1;
	double beta = alpha;
	double K1 = 0, K2 = 0;

	// RK2 ODE Solver
	for (int i = 0; i < N; i++) {
		
		// [First-point Gradient]
		K1 = odeFunc(ti, y[i]);

		// [Second-point Gradient]
		t2 = ti + alpha * h;
		y2 = y[i] + beta * K1 * h;
		K2 = odeFunc(t2, y2);

		// Estimate: y(i+1) using RK2 
		y[i + 1] = y[i] + (C1 * K1 + C2 * K2) * h;

		ti += h;

	}
}

void odeRK3(double y[], double odeFunc(const double t, const double y), const double t0, const double tf, const double h, const double y0) {

	/*Input:
	*	y[]:			A vector with the y coordinate of the solution points.
	*	odeFunc(t, y):	Name of a function file that calculates dy/dt.
	*	t0:				The first value of t.
	*	tf:				The last value of t.
	*	h:				Step size.
	*	y_init:			Initial Condition of y(1)
	*/

	// Variable Initialization
	int N = (tf - t0) / h;
	double t2 = 0, t3 = 0, t4 = 0;
	double y2 = 0, y3 = 0, y4 = 0;

	// Initial Condition
	double ti = t0;
	y[0] = y0;

	// RK Design Parameters
	double C1 = 1.0 / 6.0;
	double C2 = 4.0 / 6.0;
	double C3 = 1.0 / 6.0;

	double alpha2 = 0.5;
	double alpha3 = 1.0;
	double beta21 = 0.5;
	double beta31 = -1.0;
	double beta32 = 2.0;

	double K1 = 0;
	double K2 = 0;
	double K3 = 0;

	for (int i = 0; i < N; i++) {
		
		K1 = odeFunc(ti, y[i]);

		t2 = ti + alpha2 * h;
		y2 = y[i] + beta21 * K1 * h;
		K2 = odeFunc(t2, y2);

		t3 = ti + alpha3 * h;
		y3 = y[i] + beta31 * K1 * h + beta32 * K2 * h;
		K3 = odeFunc(t3, y3);

		y[i + 1] = y[i] + (C1 * K1 + C2 * K2 + C3 * K3) * h;

		ti += h;
	}
}

// ODE RK2 for 2nd order ODE 
void sys2RK2(double y[], double v[], void odeFuncSys(double dYdt[], const double t, const double Y[]), 
	const double t0, const double tf, const double h, const double y_init, const double v_init)
{
	/*Input:
	*	y[]:			A vector with the y coordinate of the solution points.
	*	v[]:			A vector with the v coordinate of the solution points.
	*	odeFuncSys:		Name of a function file that calculates [dy/dt ; dv/dt].
	*	t0:				The first value of t.
	*	tf:				The last value of t.
	*	h:				Step size.
	*	y_init:			The initial value of y.
	*	v_init:			The initial value of v.
	*/

	// Variable Initialization
	int N = (tf - t0) / h + 1;
	double t2 = 0;

	// Initial Condition	
	double ti = t0;
	y[0] = y_init;
	v[0] = v_init;

	// RK Design Parameters
	double alpha = 1.0;
	double beta = alpha;
	double C1 = 0.5;
	double C2 = 1.0 - C1;

	// RK Design Parameters
	double K1[2] = { 0 };           //{Ky1, Kv1}
	double K2[2] = { 0 };			//{Ky2, Kv2}
	double dYdt[2] = { 0 };
	double dYdt2[2] = { 0 };
	double Yin[2] = { 0 }; //{y2, v2}


	// RK2 ODE for 2nd ODE
	for (int i = 0; i < N - 1; i++) {

		Yin[0] = y[i];
		Yin[1] = v[i];

		odeFuncSys(dYdt, ti, Yin);
		K1[0] = dYdt[0];  //Ky1
		K1[1] = dYdt[1];  //Kv1

		t2 = ti + alpha * h;
		Yin[0] = y[i] + beta * K1[0] * h;  //y2 = y(i) + beta*Ky1*h;
		Yin[1] = v[i] + beta * K1[1] * h;  //v2 = v(i)+ beta*Kv1*h;

		odeFuncSys(dYdt2, t2, Yin);
		K2[0] = dYdt2[0];	//Ky2
		K2[1] = dYdt2[1];	//Kv2

		y[i + 1] = y[i] + (C1 * K1[0] + C2 * K2[0]) * h;
		v[i + 1] = v[i] + (C1 * K1[1] + C2 * K2[1]) * h;
		ti += h;
	}
}

void sys2RK4(double y[], double v[], void odeFuncSys(double dYdt[], const double t, const double Y[]), 
	const double t0, const double tf, const double h, const double y_init, const double v_init) {

	/*Input:
	*	y[]:			A vector with the y coordinate of the solution points.
	*	v[]:			A vector with the v coordinate of the solution points.
	*	odeFuncSys:		Name of a function file that calculates [dy/dt ; dv/dt].
	*	t0:				The first value of t.
	*	tf:				The last value of t.
	*	h:				Step size.
	*	y_init:			The initial value of y.
	*	v_init:			The initial value of v.
	*/

	// Variable Initialization
	int N = (tf - t0) / h + 1;
	double t2 = 0, t3 = 0, t4 = 0;

	// Initial Condition	
	double ti = t0;
	y[0] = y_init;
	v[0] = v_init;

	// RK Design Parameters
	double alpha2 = 0.5;
	double alpha3 = alpha2;
	double beta21 = alpha2;
	double beta32 = alpha2;

	double alpha4 = 1;
	double beta43 = alpha4;
	double beta31 = 0;
	double beta41 = beta31;
	double beta42 = beta31;

	double C1 = 1.0 / 6.0;
	double C4 = C1;
	double C2 = 2.0 / 6.0;
	double C3 = C2;

	// RK Design Parameters
	double K1[2] = { 0 };
	double K2[2] = { 0 };
	double K3[2] = { 0 };
	double K4[2] = { 0 };

	double dYdt[2] = { 0 };
	double dYdt2[2] = { 0 };
	double dYdt3[2] = { 0 };
	double dYdt4[2] = { 0 };

	double Yin[2] = { 0 };

	// RK2 ODE for 2nd ODE
	for (int i = 0; i < N - 1; i++) {

		Yin[0] = y[i];
		Yin[1] = v[i];

		//1
		odeFuncSys(dYdt, ti, Yin);
		K1[0] = dYdt[0];
		K1[1] = dYdt[1];

		//2
		t2 = ti + alpha2 * h;
		Yin[0] = y[i] + beta21 * K1[0] * h;
		Yin[1] = v[i] + beta21 * K1[1] * h;

		odeFuncSys(dYdt2, t2, Yin);
		K2[0] = dYdt2[0];
		K2[1] = dYdt2[1];

		//3
		t3 = ti + alpha3 * h;
		Yin[0] = y[i] + beta31 * K1[0] * h + beta32 * K2[0] * h;
		Yin[1] = v[i] + beta31 * K1[1] * h + beta32 * K2[1] * h;

		odeFuncSys(dYdt3, t3, Yin);
		K3[0] = dYdt3[0];
		K3[1] = dYdt3[1];

		//4
		t4 = ti + alpha4 * h;
		Yin[0] = y[i] + beta41 * K1[0] * h + beta42 * K2[0] * h + beta43 * K3[0] * h;
		Yin[1] = v[i] + beta41 * K1[1] * h + beta42 * K2[1] * h + beta43 * K3[1] * h;

		odeFuncSys(dYdt4, t4, Yin);
		K4[0] = dYdt4[0];
		K4[1] = dYdt4[1];

		//out
		y[i + 1] = y[i] + (C1 * K1[0] + C2 * K2[0] + C3 * K3[0] + C4 * K4[0]) * h;
		v[i + 1] = v[i] + (C1 * K1[1] + C2 * K2[1] + C3 * K3[1] + C4 * K4[1]) * h;
		ti += h;
	}
}

/*----------------- End ODE--------------------*/


