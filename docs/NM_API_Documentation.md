`#include "myMatrix_22000287.h"`

# SystemNonLinear

## nonlinearSys()

Non-lineaer Equation Solver

```c
Matrix nonlinearSys(Matrix Funcs(Matrix _Z), Matrix Jacob(Matrix _Z), Matrix _Z0, double tol);
```

#### **Parameters**

- **Funcs(Matrix _Z)**:
  non-linear equation that want to solve. This matrix is (n*1)
  ex)
  
  ```c
  Matrix myFuncEx1(Matrix X){
	int n = X.rows;
	Matrix F = zeros(n, 1);
	double x1 = X.at[0][0];
	double x2 = X.at[1][0];

	F.at[0][0] = x2 - 0.5 * (exp(0.5 * x1) + exp(-0.5 * x1));
	F.at[1][0] = 9 * (x1) * (x1) + 25 * x2 * x2 - 225;

	return F;}
  ```


- **Jacob(Matrix _Z)**:
  Jacobian matrix (= F'(X)). It is differential value of the matrix Funcs(Matrix _Z).
  This matrix is square (n*n)
  ex)
  ```c
  Matrix myJacobEx1(Matrix X){
	int n = X.rows;
	Matrix J = zeros(n, n);
	double x1 = X.at[0][0];
	double x2 = X.at[1][0];

	J.at[0][0] = -0.25 * (exp(0.5 * x1) - exp(-0.5 * x1));
	J.at[0][1] = 1;
	J.at[1][0] = 18 * x1;
	J.at[1][1] = 50 * x2;


	return J;}
  ```

- **_Z0**:  Matrix **_Z0** is Initial condition of Z. It is (n*1) matrix.
  ```c
  Matrix Z = zeros(n, 1);

  // Initial condition
  double z0[2] = { 2.5, 2 };
  Z = arr2Mat(z0, n, 1);
  ```

- **tol**:  Tolerance value.

#### Example code
```c
double loss = 0;
double n = 2;
Matrix J = zeros(n, n);
Matrix F = zeros(n, 1);
Matrix H = zeros(n, 1);
Matrix Z = zeros(n, 1);

// Initial condition
double z0[2] = { 2.5, 2 };
Z = arr2Mat(z0, n, 1);

Z = nonlinearSys(myFuncEx1, myJacobEx1, Z, 0.00001);
printMat(Z, "Z");
```
---

## DOF_nonlinearSys()

Non-lineaer Equation Solver for 3-DOF transformation.
angle(theta) and translation T=[dx, dy] to move the vehicle to the new pose. 
![image](https://github.com/user-attachments/assets/86119918-39c2-4fd6-84c4-9cd0dbb2dd60)

almost same as nonlinearSys().

```c
Matrix DOF_nonlinearSys(Matrix Funcs(Matrix _Z), Matrix Jacob(Matrix _Z), Matrix _Z0, double tol);
```

#### **Parameters**

- **Funcs(Matrix _Z)**:
  non-linear equation that want to solve. This matrix is (n*1)
  ex)
  
  ```c
  Matrix myFuncEx2(Matrix X){
	int n = X.rows;
	Matrix F = zeros(n, 1);
	double th = X.at[0][0];
	double dx = X.at[1][0];
	double dy = X.at[2][0];

	int P0x = 0;		int P0y = 100;
	int P1x = 0;		int P1y = -100;
	int P0x_new = 50;	double P0y_new = 186.6025;
	int P1x_new = 150;	double P1y_new = 13.3975;

	F.at[0][0] = cos(th) * P0x - sin(th) * P0y + dx - P0x_new;
	F.at[1][0] = sin(th) * P0x + cos(th) * P0y + dy - P0y_new;
	F.at[2][0] = cos(th) * P1x - sin(th) * P1y + dx - P1x_new;

	return F;}
  ```


- **Jacob(Matrix _Z)**:
  Jacobian matrix (= F'(X)). It is differential value of the matrix Funcs(Matrix _Z).
  This matrix is square (n*n)
  ex)
  ```c
  Matrix myJacobEx2(Matrix X){
	int n = X.rows;
	Matrix J = zeros(n, n);
	double th = X.at[0][0];
	double dx = X.at[1][0];
	double dy = X.at[2][0];

	int x0 = 0; 
	int y0 = 100; 
	int x1 = 0; 
	int y1 = -100;

	J.at[0][0] = (-sin(th)) * x0 - (cos(th)) * y0;
	J.at[0][1] = 1;
	J.at[0][2] = 0;

	J.at[1][0] = (cos(th)) * x0 - (sin(th)) * y0;
	J.at[1][1] = 0;
	J.at[1][2] = 1;

	J.at[2][0] = (-sin(th)) * x1 - (cos(th)) * y1;
	J.at[2][1] = 1;
	J.at[2][2] = 0;

	return J;}
  ```

- **_Z0**:  Matrix **_Z0** is Initial condition of Z. It is (n*1) matrix.
  ```c
  Matrix Z_Q2 = zeros(n, 1);

  // Initial condition
  double z0_Q2[3] = { PI * 25 / 180, 90, 90 };
  Z_Q2 = arr2Mat(z0_Q2, n_Q2, 1);
  ```

- **tol**:  Tolerance value.

#### Example code
```c
double loss_Q2 = 0;
double n_Q2 = 3;
Matrix J_Q2 = zeros(n, n);
Matrix F_Q2 = zeros(n, 1);
Matrix H_Q2 = zeros(n, 1);
Matrix Z_Q2 = zeros(n, 1);

// Initial condition
double z0_Q2[3] = { PI * 25 / 180, 90, 90 };
Z_Q2 = arr2Mat(z0_Q2, n_Q2, 1);

Z_Q2 = DOF_nonlinearSys(myFuncEx2, myJacobEx2, Z_Q2, 0.001);
Z_Q2.at[0][0] *= (180 / PI); // rad to deg
printMat(Z_Q2, "Z_Q2 ");
```
---


---


# Linear Solver

## gaussElim()

solves for vector **x** from  Ax=b,  a linear system problem  

```c
void	gaussElim(Matrix _A, Matrix _B, Matrix* _U, Matrix* _B_out);
```

#### **Parameters**

- **A**:  Matrix **A** in structure Matrix form.  Should be (nxn) square.

- **B**:  vector  **b** in structure Matrix form.  Should be (nx1) 

  

  

## solveLinear()

solves for vector **x** from  Ax=b,  a linear system problem  

```c
extern Matrix solveLinear(Matrix _A, Matrix _b, char* _method)
```

#### **Parameters**

- **A**:  Matrix **A** in structure Matrix form.  Should be (nxn) square.

- **b**:  vector  **b** in structure Matrix form.  Should be (nx1) 

- **method:  character type,** 

  - **'lu' :** LU decomposition
  - **'gauss':** Gauss elimination

  

#### Example code

```C
double A_array[] = { 1, 3, -2, 4,		2, -3, 3, -1,		-1, 7, -4, 2,		-1, 7, -4, 2 };
double b_array[] = { -11,		6,		-9,		15 };

Matrix matA = arr2Mat(A_array, M, N);
Matrix vecb = arr2Mat(b_array, M, 1);

Matrix x_lu = solveLinear(matA, vecb, "LU");
Matrix invA_gj = inv(matA, "gj");
Matrix invA_lu = inv(matA, "LU");
```





***



# Numerical Differentiation

## gradientFunc()

Solves for numerical gradient  (dy/dt) from  given equation

```c
Matrix gradientFunc(double func(double t), double xin);
```

#### Parameters

* **func**: a function as input argument 

## gradient()

Solves for numerical gradient  (dy/dt) from  a set of discrete data

```c
Matrix	gradient(Matrix _t, Matrix _y);
```

#### **Parameters**

- **t**:  vector **t** in structure Matrix form.  Should be (nx1) vector
- **y**:  vector  **y** in structure Matrix form.  Should be (nx1) vector and same length as t
- Returns **dydt** in structure Matrix form. Output is also (nx1) vector



#### Example code

```c
Matrix t = txt2Mat("", "Q1_vect");
Matrix x = txt2Mat("", "Q1_vecx");

Matrix vel = gradient(t, x);
Matrix acc = gradient(t, vel);

printMat(t, "t");
printMat(x, "x");
printMat(vel, "vel");
printMat(acc, "acc");
```

See full example code:  [TutorialDifferentiation.cpp](https://github.com/ykkimhgu/tutorial-NM/blob/main/samples/Tutorial-Differentiation.cpp)
