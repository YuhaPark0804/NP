---
description: Numerical Method example


---

# Example: API documentation

## Reference of Numerical Programming API

`#include "myNM.h"`

## Non-Linear Solver

### newtonRaphson\(\)

Solves the non-linear problem using Newton-Raphson method

```text
double newtonRaphson(double x0, double tol);
```

**Parameters**

* **x0:**  initial value.
* **tol**:  tolerance error

**Example code**

```text
double tol = 0.00001;
double x0 = 3;
double NR_result;

NR_result = newtonRaphson(x0, tol);
```

## Linear Solver

### gaussElim\(\)

solves for vector **x** from Ax=b, a linear system problem Using Gauss Elimination

```cpp
void gaussElim(Matrix _A, Matrix _B, Matrix* _U, Matrix* _B_out);
```

**Parameters**

* **A**: Matrix **A** in structure Matrix form. Should be \(nxn\) square.
* **B**: vector **b** in structure Matrix form. Should be \(nx1\)
* **U**: Matrix **U** in structure Matrix form. Should be \(nxn\) square.
* **B\_out**: vector **B\_out** in structure Matrix form. Should be \(nx1\)

**Example code**

```text
Matrix matA = txt2Mat(path, "prob1_matA");
Matrix vecb = txt2Mat(path, "prob1_vecb");
Matrix matU = zeros(matA.rows, matA.cols);
Matrix vecd = zeros(vecb.rows, vecb.cols);

gaussElim(matA, vecb, matU, vecd);
```

### inv\(\)

Find the inverse Matrix.

```text
void inv(Matrix _A, Matrix _Ainv);
```

**Parameters**

* **A**: Matrix **A** in structure Matrix form. Should be \(nxn\) square.
* **Ainv**: Matrix **Ainv** in structure Matrix form. Should be \(nxn\) square.

**Example code**

```text
Matrix matA = txt2Mat(path, "prob1_matA");
Matrix matAinv = zeros(matA.rows, matA.cols);

inv(matA, matAinv);
```

### 

## Numerical Differentiation

### gradient1D\(\)

Solve for numerical gradient \(dy/dt\) from a 1D-array form.

```text
void gradient1D(double x[], double y[], double dydx[], int m);
```

**Parameters**

* **x\[\]**: input data vector **x** in 1D-array .
* **y\[\]**: input data vector **y** in 1D-array.
* **dydx\[\]**: output vector **dydx** in 1D-array.
* **m**:  length **x** and **y**.

**Example code**

```cpp
double x[21];
    for (int i = 0; i < 21; i++) {
        x[i] = 0.2 * i;
    }
double y[] = { -5.87, -4.23, -2.55, -0.89, 0.67, 2.09, 3.31, 4.31, 5.06, 5.55, 5.78, 5.77, 5.52, 5.08, 4.46, 3.72, 2.88, 2.00, 1.10, 0.23, -0.59 };
double dydx[21];

gradient1D(x, y, dydx, 21);
```

See full example code:  [TutorialDifferentiation.cpp](https://github.com/ykkimhgu/tutorial-NM/blob/main/tutorial/Tutorial-Differentiation.cpp)





## Integration

### integral\(\)

Integral using Simpson 1/3 Method.

```text
double integral(double func(const double _x), double a, double b, int n);
```

**Parameters**

* **func**: Function **func** is defined.
* **a** is starting point of x.
* **b** is ending point of x.
* **n** is the length between **a** and **b**

**Example code**

```text
double I_simpson13 = integral(myFunc, -1, 1, 12);

double myFunc(const double _x) {
	return sqrt(1 - (_x * _x));
}
```

## ODE-IVP

### odeEU\(\)

Solve the 1st-order ODE using Euler's Explicit Method.

```text
void odeEU(double y[], double odeFunc(const double t, const double y), 
	const double t0, const double tf, const double h, const double y_init);
```

**Parameters**

- **y\[\]**: Solution of ODE in structure 1D-array form.

example code

```c
unsigned int N = ((b - a) / h) + 1;
double y_EU[200] = { 0 };	//Cannot use sytanx of y_EU[N]={0};	
```



* **odeFunc()**: Function **func** is defined.

example code

```c
// Gradient function for ODE - 1st order 
double odeFunc_rc(const double t, const double v)
{
	// Input:	 y, t
	// Output:	 dydt 

	// system modeling parameters 
	double tau = 0.01;
	double f = 100;
	double Vm = 1;
	double omega = 2 * PI * f;
	double dvdt = 0;
	
	dvdt = -v/tau + (Vm / tau) * cos(omega * t);

	return dvdt;
}
```



* **t0** is starting point.
* **tf** is ending point.
* **h** is length of step.
* **y_init** is initial value of **y\[\]**.

**Example code**

```c
double a = 0;
double b = 0.1;
double h = 0.001;

unsigned int N = ((b - a) / h) + 1;
double y_EU[200] = { 0 };

double y_init = 0;

odeEU(y_EU, odeFunc_rc, a, b, h, y_init);

double odeFunc_rc(const double t, const double v) {
	double tau = 1;
	double T = 1 / tau;
	double f = 10;
	double Vm = 1;
	double omega = 2 * PI * f;
	return  -T * v + T * Vm * cos(omega * t);
}
```

-------------------------------------------------------------------------------------------------------



### odeRK2\(\)

Solve the 1st-order ODE using Runge-Kutta 2nd order

```c
void odeRK2(double y[], double odeFunc(const double t, const double y), 
	const double t0, const double tf, const double h, const double y_init);
```

**Parameters**

- **y\[\]**: Solution of ODE in structure 1D-array form.

example code

```c
unsigned int N = ((b - a) / h) + 1;
double y_RK2[200] = { 0 };	//Cannot use sytanx of y_EU[N]={0};	
```



* **odeFunc()**: Function **func** is defined.

example code

```c
// Gradient function for ODE - 1st order 
double odeFunc_rc(const double t, const double v)
{
	// Input:	 y, t
	// Output:	 dydt 

	// system modeling parameters 
	double tau = 0.01;
	double f = 100;
	double Vm = 1;
	double omega = 2 * PI * f;
	double dvdt = 0;
	
	dvdt = -v/tau + (Vm / tau) * cos(omega * t);

	return dvdt;
}
```



* **t0** is starting point.
* **tf** is ending point.
* **h** is length of step.
* **y_init** is initial value of **y\[\]**.

**Example code**

```c
double a = 0;
double b = 0.1;
double h = 0.001;

unsigned int N = ((b - a) / h) + 1;
double y_RK2[200] = { 0 };

double y_init = 0;

//  Runge-Kutta 2nd order
odeRK2(y_RK2, odeFunc_rc, a, b, h, y_init);

double odeFunc_rc(const double t, const double v) {
	double tau = 1;
	double T = 1 / tau;
	double f = 10;
	double Vm = 1;
	double omega = 2 * PI * f;
	return  -T * v + T * Vm * cos(omega * t);
}
```

-------------------------------------------------------------------------------------------------------

## 



### odeRK3\(\)

Solve the 1st-order ODE using Runge-Kutta 3nd order

```c
void odeRK3(double y[], double odeFunc(const double t, const double y), 
	const double t0, const double tf, const double h, const double y0);
```

**Parameters**

- **y\[\]**: Solution of ODE in structure 1D-array form.

example code

```c
unsigned int N = ((b - a) / h) + 1;
double y_RK3[200] = { 0 };	//Cannot use sytanx of y_EU[N]={0};	
```



* **odeFunc()**: Function **func** is defined.

example code

```c
// Gradient function for ODE - 1st order 
double odeFunc_rc(const double t, const double v)
{
	// Input:	 y, t
	// Output:	 dydt 

	// system modeling parameters 
	double tau = 0.01;
	double f = 100;
	double Vm = 1;
	double omega = 2 * PI * f;
	double dvdt = 0;
	
	dvdt = -v/tau + (Vm / tau) * cos(omega * t);

	return dvdt;
}
```



* **t0** is starting point.
* **tf** is ending point.
* **h** is length of step.
* **y0** is initial value of **y\[\]**.

**Example code**

```c
double a = 0;
double b = 0.1;
double h = 0.001;

unsigned int N = ((b - a) / h) + 1;
double y_RK3[200] = { 0 };

double y_init = 0;

//  Runge-Kutta 3rd order
odeRK3(y_RK3, odeFunc_rc, a, b, h, y_init);

printf("----------------------------------------------------------------\n");
printf("			       1st ODE - IVP  Results						\n");
printf("----------------------------------------------------------------\n");

printf("i\t t\t\t yEU\t\t yRK2\t\t yRK3\t \n\n");
for (int i = 0; i < N; i++)
	printf("%d\t %0.5f\t %.5f\t %.5f\t %.5f\t \n", i, a + i * h, y_EU[i], y_RK2[i], y_RK3[i]);
printf("\n");

double odeFunc_rc(const double t, const double v) {
	double tau = 1;
	double T = 1 / tau;
	double f = 10;
	double Vm = 1;
	double omega = 2 * PI * f;
	return  -T * v + T * Vm * cos(omega * t);
}
```

-------------------------------------------------------------------------------------------------------

## 



### sys2RK2\(\)

Solve the 2st-order ODE using Runge-Kutta 2nd-order

```c
void sys2RK2(double y[], double v[], void odeFuncSys(double dYdt[], const double t, const double Y[]), 
	const double t0, const double tf, const double h, const double y_init, const double v_init);
```

**Parameters**

- **y\[\]**: Solution of ODE in structure 1D-array form. (y' = z(t))
- **v\[\]**: Solution of ODE in structure 1D-array form. (y'' = z')

* **odeFuncSys()**: Function **func** is defined.

example code

```c
// Gradient function for ODE - 2nd order 
void odeFunc_mck(double dYdt[], const double t, const double Y[])
{
	// Input:	 vecY = [y; v], t
	// Output:	 dYdt = [dydt; dvdt]

	// system modeling parameters 
	double	m		= 1.0;					//[kg]
	double	k		= 6.9;					//[N/m]
	double	c		= 7.0;				    //[N/m/s]
	double	A		= 2.0;					//[N]
	double	f		= 5.0;					//[Hz]
	double	omega	= 2.0 * PI * f;
	double	Fin		= A * cos(omega * t);

	// output
	dYdt[0] = Y[1];  //dydt = v
	dYdt[1] = (1.0 / m) * (-k * Y[0] - c * Y[1] + Fin);
}
```



* **t0** is starting point.
* **tf** is ending point.
* **h** is length of step.
* **y_init** is initial value of **y\[\]**.
* **z_init** is initial value of **z\[\]**.

**Example code**

```c
//Parameter Definitions
double t0 = 0;
double tf = 1;
h = 0.01;
N = (tf - t0) / h + 1;

// Runge-Kutta 2nd order
double Y_2RK2[200] = { 0 };
double V_2RK2[200] = { 0 };

// Initial values
double y_initQ2 = 0;
double v_initQ2 = 0.2;

//  Runge-Kutta 2nd order
sys2RK2(Y_2RK2, V_2RK2, odeFunc_mck, t0, tf, h, y_initQ2, v_initQ2);

// Gradient function for ODE - 2nd order 
void odeFunc_mck(double dYdt[], const double t, const double Y[])
{
	// Input:	 vecY = [y; v], t
	// Output:	 dYdt = [dydt; dvdt]

	// system modeling parameters 
	double	m		= 1.0;					//[kg]
	double	k		= 6.9;					//[N/m]
	double	c		= 7.0;				    //[N/m/s]
	double	A		= 2.0;					//[N]
	double	f		= 5.0;					//[Hz]
	double	omega	= 2.0 * PI * f;
	double	Fin		= A * cos(omega * t);

	// output
	dYdt[0] = Y[1];  //dydt = v
	dYdt[1] = (1.0 / m) * (-k * Y[0] - c * Y[1] + Fin);
}
```

-------------------------------------------------------------------------------------------------------







### sys2RK4\(\)

Solve the 2st-order ODE using Runge-Kutta 4nd-order

```c
void sys2RK4(double y[], double v[], void odeFuncSys(double dYdt[], const double t, const double Y[]), 
	const double t0, const double tf, const double h, const double y_init, const double v_init);
```

**Parameters**

- **y\[\]**: Solution of ODE in structure 1D-array form. (y' = z(t))
- **v\[\]**: Solution of ODE in structure 1D-array form. (y'' = z')

* **odeFuncSys()**: Function **func** is defined.

example code

```c
// Gradient function for ODE - 2nd order 
void odeFunc_mck(double dYdt[], const double t, const double Y[])
{
	// Input:	 vecY = [y; v], t
	// Output:	 dYdt = [dydt; dvdt]

	// system modeling parameters 
	double	m		= 1.0;					//[kg]
	double	k		= 6.9;					//[N/m]
	double	c		= 7.0;				    //[N/m/s]
	double	A		= 2.0;					//[N]
	double	f		= 5.0;					//[Hz]
	double	omega	= 2.0 * PI * f;
	double	Fin		= A * cos(omega * t);

	// output
	dYdt[0] = Y[1];  //dydt = v
	dYdt[1] = (1.0 / m) * (-k * Y[0] - c * Y[1] + Fin);
}
```



* **t0** is starting point.
* **tf** is ending point.
* **h** is length of step.
* **y_init** is initial value of **y\[\]**.
* **z_init** is initial value of **z\[\]**.

**Example code**

```c
//Parameter Definitions
double t0 = 0;
double tf = 1;
h = 0.01;
N = (tf - t0) / h + 1;

// Runge-Kutta 4th order
double Y_2RK4[200] = { 0 };
double V_2RK4[200] = { 0 };

// Initial values
double y_initQ2 = 0;
double v_initQ2 = 0.2;

// Runge-Kutta 4th order
sys2RK4(Y_2RK4, V_2RK4, odeFunc_mck, t0, tf, h, y_initQ2, v_initQ2);

printf("----------------------------------------------------------------\n");
printf("			       2nd order ODE - IVP  Results					\n");
printf("----------------------------------------------------------------\n");

printf("i\t t\t\t Y_2RK2\t\t V_2RK2\t\t Y_2RK4\t\t V_2RK4\t\t \n\n");
for (int i = 0; i < N; i++)
	printf("%d\t %0.5f\t %.5f\t %.5f\t %.5f\t %.5f\t\n", i, a + i * h, Y_2RK2[i], V_2RK2[i], Y_2RK4[i], V_2RK4[i]);
printf("\n");


// Gradient function for ODE - 2nd order 
void odeFunc_mck(double dYdt[], const double t, const double Y[])
{
	// Input:	 vecY = [y; v], t
	// Output:	 dYdt = [dydt; dvdt]

	// system modeling parameters 
	double	m		= 1.0;					//[kg]
	double	k		= 6.9;					//[N/m]
	double	c		= 7.0;				    //[N/m/s]
	double	A		= 2.0;					//[N]
	double	f		= 5.0;					//[Hz]
	double	omega	= 2.0 * PI * f;
	double	Fin		= A * cos(omega * t);

	// output
	dYdt[0] = Y[1];  //dydt = v
	dYdt[1] = (1.0 / m) * (-k * Y[0] - c * Y[1] + Fin);
}
```

-------------------------------------------------------------------------------------------------------





### Pred_Corr\(\)

Solve ODE using explicit(predictor) and implicit(corrector)

predictor: Explicit formula to estimate the solution `y_i+1`

Corrector: Use `y_i+1`(predictor) on implicit formula to obtain more accurate `y_i+1`

```c
void Pred_Corr(double y[], double odeFunc(const double t, const double y),
	const double t0, const double tf, const double h, const double y_init);
```

**Parameters**

- **y\[\]**: Solution of ODE in structure 1D-array form. (y' = z(t))
- **v\[\]**: Solution of ODE in structure 1D-array form. (y'' = z')

* **odeFuncSys()**: Function **func** is defined.

* **t0** is starting point.
* **tf** is ending point.
* **h** is length of step.
* **y_init** is initial value of **y\[\]**.

**Example code**

```c
// Initial Conditions
double a = 0;
double b = 0.1;
double h = 0.001;

unsigned int N = ((b - a) / h) + 1;

double y_PC[200]  = { 0 };

// IVP initial values
double y_init = 0;

//Predictor-Corrector methods
Pred_Corr(y_PC, odeFunc_rc, a, b, h, y_init);

// Print outputs
printf("----------------------------------------------------------------\n");
printf("			       1st ODE - IVP  Results						\n");
printf("----------------------------------------------------------------\n");

printf("i\t t\t\t yEU\t\t yRK2\t\t yRK3\t\t yPC\t\t yAB\t \n\n");
for (int i = 0; i < N; i++)
	printf("%d\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t \n", i, a + i * h, y_EU[i], y_RK2[i], y_RK3[i], y_PC[i], y_AB[i]);
printf("\n");


// Gradient function for ODE - 1st order 
double odeFunc_rc(const double t, const double v)
{
	// Input:	 y, t
	// Output:	 dydt 

	// system modeling parameters 
	double tau = 0.01;
	double f = 100;
	double Vm = 1;
	double omega = 2 * PI * f;
	double dvdt = 0;
	
	dvdt = -v/tau + (Vm / tau) * cos(omega * t);

	return dvdt;
}
```

-------------------------------------------------------------------------------------------------------



### Adam_Bash3\(\)

Solve 1st order ODE using explicit multi step

```c
void Adam_Bash3(double y[], double odeFunc(const double t, const double y), 
	const double t0, const double tf, const double h, const double y_init);
```

**Parameters**

- **y\[\]**: Solution of ODE in structure 1D-array form. (y' = z(t))
- **v\[\]**: Solution of ODE in structure 1D-array form. (y'' = z')

* **odeFuncSys()**: Function **func** is defined.

* **t0** is starting point.
* **tf** is ending point.
* **h** is length of step.
* **y_init** is initial value of **y\[\]**.

**Example code**

```c
// Initial Conditions
double a = 0;
double b = 0.1;
double h = 0.001;

unsigned int N = ((b - a) / h) + 1;

double y_AB[200] = { 0 };

// IVP initial values
double y_init = 0;

//Adam-Bashforth 3-order methods
Adam_Bash3(y_AB, odeFunc_rc, a, b, h, y_init);

// Print outputs
printf("----------------------------------------------------------------\n");
printf("			       1st ODE - IVP  Results						\n");
printf("----------------------------------------------------------------\n");

printf("i\t t\t\t yEU\t\t yRK2\t\t yRK3\t\t yPC\t\t yAB\t \n\n");
for (int i = 0; i < N; i++)
	printf("%d\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t \n", i, a + i * h, y_EU[i], y_RK2[i], y_RK3[i], y_PC[i], y_AB[i]);
printf("\n");


// Gradient function for ODE - 1st order 
double odeFunc_rc(const double t, const double v)
{
	// Input:	 y, t
	// Output:	 dydt 

	// system modeling parameters 
	double tau = 0.01;
	double f = 100;
	double Vm = 1;
	double omega = 2 * PI * f;
	double dvdt = 0;
	
	dvdt = -v/tau + (Vm / tau) * cos(omega * t);

	return dvdt;
}
```







