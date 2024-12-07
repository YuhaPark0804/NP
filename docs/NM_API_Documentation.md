`#include "myMatrix_22000287.h"`

## Eigen Value

### MatrixNorm\(\)

Solve to find the norm  of Matrix.

Matrix A = [a0, a1;...];

||n|| = sqrt(a0^2 + a1^2 + ... + an^2)

```c
double MatrixNorm(Matrix& A);
```

**Parameters**

- **A**: Matrix that wants to find the norm.



**Example code**

```c
Matrix c = createMat(n, 1);
MatrixNorm(c);
```



-------------------------------------------------------------------------------------------------------



### QRdecomp\(\)

To solve Eigenvalue & Eigenvector, A has to be decomposition (Matrix A = QR) 

Using HouseHold Matrix

Transform A -> U, while preserving eigenvalues. (use similar matrix)

```c
void QRdecomp(Matrix& A, Matrix& Q, Matrix& R);
```

**Parameters**

- **A**: Matrix that wants to decompose to A=QR. A: (n*n)
- **Q**: Matrix that orthonormal. (inv(Q) = transpose(Q)) Q: (n*n)
- **R**: Matrix that Upper triangular. R (n*n)



**Example code**

```c
Matrix U = copyMat(A);
Matrix Q = eye(m, n);
Matrix R = copyMat(A);
QRdecomp(U, Q, R); // update Q, R
```



-------------------------------------------------------------------------------------------------------





### eigval\(\)

To solve Eigenvalue. (Estimating Eigenvalues using QR factorization and Iteration)

Find a similar matrix U that preserves the eigenvalues.

```c
Matrix eigval(Matrix& A);
```

**Parameters**

- **A**: Matrix that wants to find eigenvalue.

**Output**

- **lamda**: (n*1) matrix



**Example code**

```c
Matrix matA_Q1 = txt2Mat(path, "prob1_matA");
Matrix eigVals_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);

eigVals_Q1 = eigval(matA_Q1);

printf("\n[Eigen value]\n\n");
printMat(eigVals_Q1,"eigVals_Q1");
```



-------------------------------------------------------------------------------------------------------







### eigvec\(\)

To solve Eigen vector.

```c
Matrix eigvec(Matrix A);
```

**Parameters**

- **A**: Matrix that wants to find eigen vector. (n*n) Matrix



**Output**

- **V**: (n*n) Matrix



**Example code**

```c
Matrix matA_Q1 = txt2Mat(path, "prob1_matA");
Matrix eigVals_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);

eigVecs_Q1 = eigvec(matA_Q1);

printf("\n[Eigen vector]\n\n");
printMat(eigVecs_Q1,"eigVecs_Q1");
```



-------------------------------------------------------------------------------------------------------





### eig\(\)

To solve Eigen vector & Eigen value.

```c
void eig(Matrix A, Matrix& V, Matrix& D);
```

**Parameters**

- **A**: Matrix that wants to find the eigenvector and eigenvalue.
- **V**: Matrix that wants to find eigenvector.
- **D**: Matrix that wants to find eigenvalue.



**No Output, Just update**

- **V**: eigenvector. (n*n) Matrix
- **D**: eigenvalue. (n*n) Matrix; eye Matrix but pivot is eigenvalue. 



**Example code**

```c
Matrix matA_Q1 = txt2Mat(path, "prob1_matA");
Matrix matD_Q1 = eye(matA_Q1.rows, matA_Q1.cols);
Matrix matV_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);

eig(matA_Q1, matV_Q1, matD_Q1);

printf("\n[Eigen value]\n\n");
printMat(matD_Q1, "matD_Q1");

printf("\n[Eigen vector]\n\n");
printMat(matV_Q1, "matV_Q1");
```



-------------------------------------------------------------------------------------------------------














## arr2Mat()

// Create a Matrix from 1D-array
if arr2Mat(1Darray, rows = 2, cols = 2)
	[data1, data2, data3, data4] (4*1) -> [data1, data2; data3, data4] (2*2)

```c
Matrix	arr2Mat(double* _1Darray, int _rows, int _cols);
```

#### **Parameters**

- **_1Darray**:
  1D array that want to change 

- **_rows**:
  rows that want to make.

- **_cols**:
  cols that want to make.

#### Example code

```c
double T[] = { 30, 40, 50, 60, 70, 80 };
double P[] = { 1.05, 1.07, 1.09, 1.14, 1.17, 1.21 };
int m_Q1 = 6;	// length of dataset

Matrix matT = arr2Mat(T, m_Q1, 1);
Matrix matP = arr2Mat(P, m_Q1, 1);
```

---

# curvefit

## linearFit_mat()

LinearRegration calculates the coefficients a1 and a0 of the linear
equation y = a1*x + a0 that best fit n = 1 data points.
// Calculates coefficients of least squares regression - Line

```c
Matrix	linearFit_mat(Matrix _X, Matrix _Y);
```

#### **Parameters**

- **_X**:
  A vector with the coordinates x of the data points.

- **_Y**:
  A vector with the coordinates y of the data points.

- Output variable:
  - a1   The coefficient a1.
  - a0   The coefficient a0.

#### Example code

```c
// Initial Conditions
double T[] = { 30, 40, 50, 60, 70, 80 };
double P[] = { 1.05, 1.07, 1.09, 1.14, 1.17, 1.21 };
double Z_Q1[2] = { 0 };
int orderN = 1;		// nth order
int m_Q1 = 6;		// length of dataset
int mx = 0; int my = 0;


// Check Is length(X)~= length(Y) ? Exit: Continue 
mx = sizeof(T) / sizeof(T[0]);
my = sizeof(P) / sizeof(P[0]);
if (my != mx) {
	printf("[ERROR] X and Y have different length!! \n\r");
	return 0;
}

Matrix matT = arr2Mat(T, m_Q1, 1);
Matrix matP = arr2Mat(P, m_Q1, 1);
Matrix vecZ = linearFit_mat(matT, matP);

printf("Z_Q1= [a0,  a1] \n\r");
printMat(vecZ_Q1, "Z_Q1");
```

---

## polyFit_mat()

LinearRegration calculates the coefficients a0 to an of the polynomial
equation y = a0 + a1*x + ... + an*(x^n) that best fit n data points.
// Calculates coefficients of least squares regression - Nth order polynomial

```c
Matrix	polyFit_mat(Matrix _vecX, Matrix _vecY, int orderN);
```

#### **Parameters**

- **_vecX**:
  A vector with the coordinates x of the data points.

- **_vecY**:
  A vector with the coordinates y of the data points.

- **orderN**:
  Order of polynomial

- Output variable:
  - Pcoef = [an ... a4 a3 a2 a1 a0] 


#### Example code

```c
// Initial Conditions
int m_Q2 = 16;			// data length
double Stress[] = { 0, 3, 4.5, 5.8, 5.9, 5.8, 6.2, 7.4, 9.6, 15.6, 20.7, 26.7,31.1, 35.6, 39.3, 41.5 };
double Strain[16] = { 0 };
for (int k = 0; k < m_Q2; k++)
	Strain[k] = 0.4 * k;


// Is length(X)~= length(Y) ? Exit: Continue 
mx = sizeof(Stress) / sizeof(Stress[0]);
my = sizeof(Strain) / sizeof(Strain[0]);

printf("Q2: mx=%d \t my=%d \n\r ", mx, my);
if (my != mx) {
	printf("[ERROR] X and Y have different length!! \n\r");
	return 0;
}

orderN = 4;	// nth order

Matrix matStrain = arr2Mat(Strain, m_Q2, 1);
Matrix matStress = arr2Mat(Stress, m_Q2, 1);
Matrix vecZ_Q2 = polyFit_mat(matStrain, matStress, orderN);

printf("----------------------------------------------------------------\n");
printf("\t\t Part 1-2: Polynomial  Fitting								\n");
printf("----------------------------------------------------------------\n");
printf("Z_Q2= [a0,  a1, .. a4] \n\r");
printMat(vecZ_Q2, "Z_Q2");
```

---

## expFit_mat()

LinearRegration calculates the coefficients a1 and a0 of the linear
equation y = b*exp^(mx) that best fit n data points.
Linear form: ln(y) = m*x + ln(b)
// Calculates coefficients of least squares regression - Exponential Curve-fit

```c
Matrix	expFit_mat(Matrix _X, Matrix _Y);
```

#### **Parameters**

- **_X**:
  A vector with the coordinates x of the data points.

- **_Y**:
  A vector with the coordinates y of the data points.

- Output variable:
  - c_hat = [c0, c1, ..., cn] 

![image](https://github.com/user-attachments/assets/a51f2ca7-71f0-4bbc-b78f-232bca0ac743)

#### Example code

```c
int m_Q3 = 15;				// data length
double Voltage[] = { 9.7, 8.1, 6.6, 5.1, 4.4, 3.7, 2.8, 2.4, 2.0, 1.6, 1.4, 1.1, 0.85, 0.69, 0.6 };
double Time[15] = { 0 };

double Sum = 0;
for (int k = 0; k < m_Q3; k++) {
	Time[k] = 2.0 * (k + 1);
}

	
Matrix matTime    = arr2Mat(Time, m_Q3, 1);
Matrix matVoltage = arr2Mat(Voltage, m_Q3, 1);
Matrix c_hat_Q3   = expFit_mat(matTime, matVoltage);

printf("----------------------------------------------------------------\n");
printf("\t\t Part 2: Exponential Fitting								\n");
printf("----------------------------------------------------------------\n");
printf("c_hat_Q3 = [c0,  c1, .. c4] \n\r");
printMat(c_hat_Q3, "c_hat_Q3");

/*==========================================================================*/
/*					Apply your numerical method algorithm					*/
/*==========================================================================*/

	
double c0 = c_hat_Q3.at[0][0];
double c1 = c_hat_Q3.at[1][0];
double R  = 5e+06;
double V0 = exp(c0);
double C  = -1 / (R * c1);

printf("V0: %.4f, C: %.4fe-06\n", V0, C*1e+06);
```

---



---

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
  ```


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
void gaussElim(Matrix A, Matrix b, Matrix& U, Matrix& d);
```

#### **Parameters**

- **A**:  Matrix **A** in structure Matrix form.  Should be (nxn) square.

- **b**:  vector  **b** in structure Matrix form.  Should be (nx1) 

- **U**:  Matrix  **U** in structure Matrix form. Upper triangular matrix. Should be (nxn) 

- **d**:  vector  **d** in structure Matrix form.  Should be (nx1) 

  **Example code**

  ```c
  Matrix matA_Q1 = txt2Mat(path, "prob1_matA");
  Matrix vecb_Q1 = txt2Mat(path, "prob1_vecb");
  
  Matrix matU_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);
  Matrix matd_Q1 = createMat(vecb_Q1.rows, vecb_Q1.cols);
  Matrix matx_Q1 = createMat(vecb_Q1.rows, vecb_Q1.cols);
  
  printMat(matA_Q1, "matA_Q1");
  printMat(vecb_Q1, "vecb_Q1");
  
  gaussElim(matA_Q1, vecb_Q1, matU_Q1, matd_Q1); //use this
  backsub(matU_Q1, matd_Q1, matx_Q1);
  
  printMat(matU_Q1, "matU_Q1");
  printMat(matd_Q1, "matd_Q1");
  printMat(matx_Q1, "matx_Q1");
  ```

  

  

  ---

  ## backsub()

  solves for vector **x** from  Ux=y,  a linear system problem . 

  ```c
  void backsub(Matrix U, Matrix& y, Matrix& x);
  ```

  #### **Parameters**

  - **U**:  Matrix **U** in structure Matrix form.  Should be (nxn) square.

  - **y**:  vector  **y** in structure Matrix form.  Should be (nx1) 

  - **x**:  vector  **x** in structure Matrix form. Should be (nx1) 

    

    **Example code**

    ```c
    Matrix matA_Q1 = txt2Mat(path, "prob1_matA");
    Matrix vecb_Q1 = txt2Mat(path, "prob1_vecb");
    
    Matrix matU_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);
    Matrix matd_Q1 = createMat(vecb_Q1.rows, vecb_Q1.cols);
    Matrix matx_Q1 = createMat(vecb_Q1.rows, vecb_Q1.cols);
    
    printMat(matA_Q1, "matA_Q1");
    printMat(vecb_Q1, "vecb_Q1");
    
    gaussElim(matA_Q1, vecb_Q1, matU_Q1, matd_Q1);
    backsub(matU_Q1, matd_Q1, matx_Q1);
    
    printMat(matU_Q1, "matU_Q1");
    printMat(matd_Q1, "matd_Q1");
    printMat(matx_Q1, "matx_Q1");
    ```

    





---

## fwdsub()

solves for vector **y** from  Ly=b,  a linear system problem . 

```c
void fwdsub(Matrix& L, Matrix& b, Matrix& y);
```

#### **Parameters**

- **L**:  Matrix **L** in structure Matrix form.  Should be (nxn) square.

- **y**:  vector  **y** in structure Matrix form.  Should be (nx1) 

- **b**:  vector  **b** in structure Matrix form. Should be (nx1) 

  

  **Example code**

  ```c
  void solveLU(Matrix& L, Matrix& U, Matrix& P, Matrix& b, Matrix& x) {  
  
  	//PAx=Pb
  	//LUx=Pb
  	//Ly=b*
  	//Ux=y
  	
  	Matrix y = createMat(x.rows, x.cols);
  	Matrix Pb = multMat(P, b);
  
  	fwdsub(L, Pb, y);  //Ly=b*
  	backsub(U, y, x); //Ux=y
  		
  	freeMat(y);
  	freeMat(Pb);
  }
  ```

  







---

## gaussElim_pivot()

solves LU decomposition of A with pivoting in terms of PA=LU

```c
void gaussElim_pivot(Matrix A, Matrix b, Matrix& U, Matrix& P, Matrix& L, Matrix& d);
```

#### **Parameters**

- **A**:  Matrix **A** in structure Matrix form.  Should be (nxn) square.

- **b**:  vector  **b** in structure Matrix form.  Should be (nx1) 

- **U**:  Matrix  **U** in structure Matrix form. CopyMatrix of A. Should be (nxn) 

- **P**:  Matrix  **P** in structure Matrix form. eye matrix. Should be (nxn) 

- **L**:  Matrix  **L** in structure Matrix form. eye matrix. Should be (nxn) 

- **d**:  vector  **d** in structure Matrix form.  CopyMatrix of b. Should be (nx1) 

  

  **Example code**

  ```c
  // Option 1:  Read from datafile
  Matrix matA_Q1 = txt2Mat(path, "prob1_matA");
  Matrix vecb_Q1 = txt2Mat(path, "prob1_vecb");
  
  // Option 2:  Create a zero matrix with specific size
  Matrix matU_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);
  Matrix matd_Q1 = createMat(vecb_Q1.rows, vecb_Q1.cols);
  Matrix matx_Q1 = createMat(vecb_Q1.rows, vecb_Q1.cols);
  
  Matrix matP_test = eye(test_matA.rows, test_matA.cols); 
  Matrix matL_test = createMat(test_matA.rows, test_matA.cols); 
  
  printMat(matA_Q1, "matA_Q1");
  printMat(vecb_Q1, "vecb_Q1");
  
  //-----------------//
  gaussElim_pivot(matA_Q1, vecb_Q1, matU_Q1, matP_test, matL_test, matd_Q1);
  //-----------------//
  backsub(matU_Q1, matd_Q1, matx_Q1);
  
  printMat(matU_Q1, "matU_Q1");
  printMat(matd_Q1, "matd_Q1");
  printMat(matx_Q1, "matx_Q1");
  ```







---

## LUdecomp()

solves LU decomposition of A with pivoting in terms of PA=LU.

pivoting: divided by max term 

```c
void LUdecomp(Matrix A, Matrix& L, Matrix& U, Matrix& P);
```

#### **Parameters**

- **A**:  Matrix **A** in structure Matrix form.  Should be (nxn) square.

- **L**:  Matrix  **L** in structure Matrix form. eye matrix. Should be (nxn) 

- **U**:  Matrix  **U** in structure Matrix form. CopyMatrix of A. Should be (nxn) 

- **P**:  Matrix  **P** in structure Matrix form. eye matrix. Should be (nxn) 

  

  **Example code**

  ```c
  // Option 1:  Read from datafile
  Matrix test_matA = txt2Mat(path, "test_matA"); 
  
  // Option 2:  Create a zero matrix with specific size
  Matrix matU_test = createMat(test_matA.rows, test_matA.cols);
  Matrix matP_test = eye(test_matA.rows, test_matA.cols); 
  Matrix matL_test = createMat(test_matA.rows, test_matA.cols); 
  Matrix matx_test = createMat(vecb_Q1.rows, vecb_Q1.cols); 
  Matrix maty_test = createMat(vecb_Q1.rows, vecb_Q1.cols); 
  
  
  printMat(test_matA, "test_matA"); 
  
  //-----------------//
  LUdecomp(test_matA, matL_test, matU_test, matP_test);
  //-----------------//
  
  printMat(matU_test, "matU_test");
  printMat(matL_test, "matL_test");
  printMat(matP_test, "matP_test");
  ```



**Example code 2**

```c
// Option 1:  Read from datafile
Matrix matA_Q1 = txt2Mat(path, "prob1_matA");
Matrix vecb_Q1 = txt2Mat(path, "prob1_vecb"); 

// Option 2:  Create a zero matrix with specific size
Matrix matU_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);
Matrix matd_Q1 = createMat(vecb_Q1.rows, vecb_Q1.cols);
Matrix matx_Q1 = createMat(vecb_Q1.rows, vecb_Q1.cols);
Matrix matL_Q1 = zeros(matA_Q1.rows, matA_Q1.cols); 
Matrix matP_Q1 = eye(matA_Q1.rows, matA_Q1.cols);
Matrix matAinv_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);

/// Q1
printMat(matA_Q1, "matA_Q1");
printMat(vecb_Q1, "vecb_Q1");

LUdecomp(matA_Q1, matL_Q1, matU_Q1, matP_Q1);
solveLU(matL_Q1, matU_Q1, matP_Q1, vecb_Q1, matx_Q1);

//invMat(matA_Q1, matAinv_Q1);            //Find inverse mat A
//matx_Q1 = multMat(matAinv_Q1, vecb_Q1); //x= (A^-1) * b

printMat(matU_Q1, "matU_Q1");
printMat(matx_Q1, "matx_Q1");
```



---







---

## solveLU()

After LU decomposition of A, 

do forward sub + backward sub.

```c
void solveLU(Matrix& L, Matrix& U, Matrix& P, Matrix& b, Matrix& x);
```

#### **Parameters**

- **L**:  Matrix  **L** in structure Matrix form. Should be (nxn) 
- **U**:  Matrix  **U** in structure Matrix form. Should be (nxn) 
- **P**:  Matrix  **P** in structure Matrix form. eye matrix. Should be (nxn) 

- **b**:  vector  **b** in structure Matrix form.  Should be (nx1) 

- **x**:  vector  **x** in structure Matrix form. Should be (nx1) 

  

  **Example code**

  ```c
  // Option 1:  Read from datafile
  Matrix matA_Q1 = txt2Mat(path, "prob1_matA");
  Matrix vecb_Q1 = txt2Mat(path, "prob1_vecb"); 
  
  // Option 2:  Create a zero matrix with specific size
  Matrix matU_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);
  Matrix matd_Q1 = createMat(vecb_Q1.rows, vecb_Q1.cols);
  Matrix matx_Q1 = createMat(vecb_Q1.rows, vecb_Q1.cols);
  Matrix matL_Q1 = zeros(matA_Q1.rows, matA_Q1.cols); 
  Matrix matP_Q1 = eye(matA_Q1.rows, matA_Q1.cols);
  Matrix matAinv_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);
  
  /// Q1
  printMat(matA_Q1, "matA_Q1");
  printMat(vecb_Q1, "vecb_Q1");
  
  LUdecomp(matA_Q1, matL_Q1, matU_Q1, matP_Q1);
  solveLU(matL_Q1, matU_Q1, matP_Q1, vecb_Q1, matx_Q1);//////
  
  //invMat(matA_Q1, matAinv_Q1);            //Find inverse mat A
  //matx_Q1 = multMat(matAinv_Q1, vecb_Q1); //x= (A^-1) * b
  
  printMat(matU_Q1, "matU_Q1");
  printMat(matx_Q1, "matx_Q1");
  ```









---

## invMat()

Make inverse Matrix.

```c
double invMat(Matrix A, Matrix& Ainv);
```

#### **Parameters**

- **A**:  Matrix  **A** in structure Matrix form. Should be (nxn) 

- **Ainv**:  Matrix  **Ainv** in structure Matrix form. Should be (nxn) 

  

**Example code**

```c
// Option 1:  Read from datafile
Matrix matA_Q1 = txt2Mat(path, "prob1_matA");
Matrix vecb_Q1 = txt2Mat(path, "prob1_vecb"); 

// Option 2:  Create a zero matrix with specific size
Matrix matU_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);
Matrix matd_Q1 = createMat(vecb_Q1.rows, vecb_Q1.cols);
Matrix matx_Q1 = createMat(vecb_Q1.rows, vecb_Q1.cols);
Matrix matL_Q1 = zeros(matA_Q1.rows, matA_Q1.cols); 
Matrix matP_Q1 = eye(matA_Q1.rows, matA_Q1.cols);
Matrix matAinv_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);

/// Q1
printMat(matA_Q1, "matA_Q1");
printMat(vecb_Q1, "vecb_Q1");

//LUdecomp(matA_Q1, matL_Q1, matU_Q1, matP_Q1);
//solveLU(matL_Q1, matU_Q1, matP_Q1, vecb_Q1, matx_Q1);//////

invMat(matA_Q1, matAinv_Q1);            //Find inverse mat A
matx_Q1 = multMat(matAinv_Q1, vecb_Q1); //x= (A^-1) * b

//printMat(matU_Q1, "matU_Q1");
printMat(matx_Q1, "matx_Q1");
```



---

# fin



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
