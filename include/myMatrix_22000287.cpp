/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : Yuha Park
Created          : 15-10-2024
Modified         : 27-10-2024
Language/ver     : C++ in MSVS2019

Description      : myMatrix.cpp
----------------------------------------------------------------*/

#include "myMatrix_22000287.h"



// Free a memory allocated matrix
void	freeMat(Matrix _A)
{
	// 1. Free allocated column memory
	for (int i = 0; i < _A.rows; i++)
		free(_A.at[i]);
	// 2. Free allocated row memory
	free(_A.at);
}

// Create a matrix from a text file
Matrix	txt2Mat(std::string _filePath, std::string _fileName)
{
	std::ifstream file;
	std::string temp_string, objFile = _filePath + _fileName + ".txt";
	int temp_int = 0, nRows = 0;

	file.open(objFile);
	if (!file.is_open()) {
		printf("\n*********************************************");
		printf("\n  Could not access file: 'txt2Mat' function");
		printf("\n*********************************************\n");
		return createMat(0, 0);
	}
	while (getline(file, temp_string, '\t'))
		temp_int++;
	file.close();

	file.open(objFile);
	while (getline(file, temp_string, '\n'))
		nRows++;
	file.close();

	int nCols = (temp_int - 1) / nRows + 1;
	Matrix Out = createMat(nRows, nCols);

	file.open(objFile);
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < nCols; j++) {
			file >> temp_string;
			Out.at[i][j] = stof(temp_string);
		}
	file.close();

	return Out;
}

// Create Matrix with specified size
Matrix	createMat(int _rows, int _cols)
{
	// check matrix dimension
	if (_rows < 0 || _cols < 0) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'createMat' function");
		printf("\n****************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out;
	// 1. Allocate row array first
	Out.at = (double**)malloc(sizeof(double*) * _rows);
	// 2. Then, allocate column 
	for (int i = 0; i < _rows; i++)
		Out.at[i] = (double*)malloc(sizeof(double) * _cols);
	// 3. Initialize row & column values of a matrix
	Out.rows = _rows;
	Out.cols = _cols;

	// 4. Initialize with zero (optional)
	initMat(Out, 0);
	return Out;
}


// initialization of Matrix elements
void	initMat(Matrix _A, double _val)
{
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			_A.at[i][j] = _val;
}

// Print matrix
void	printMat(Matrix _A, const char* _name)
{
	printf("%s =\n", _name);
	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++)
			printf("%15.4f\t", _A.at[i][j]);
		printf("\n");
	}
	printf("\n");
}

// Matrix addition
Matrix	addMat(Matrix _A, Matrix _B)
{
	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_A.rows, _B.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j] = _A.at[i][j] + _B.at[i][j];

	return Out;
}






//////////////////////////////////////////////////////////////////
/*				Tutorial	&  Assignment						*/
//////////////////////////////////////////////////////////////////

// Create matrix of all zeros
extern Matrix	zeros(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);

	initMat(Out, 0);

	return Out;
}


// Create matrix of all ones
extern Matrix	ones(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);
	
	initMat(Out, 1);

	return Out;
}


// Create identity matrix
// Assume square matrix (rows == cols)
extern Matrix eye(int _rows, int _cols)
{

	// error HAndling
	// cheack if (rows == cols)
	// Give warning messag
	// Exit.

	if (_rows != _cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: NOT a Square Matrix (rows != cols)");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_rows, _cols);
	
	for (int i = 0; i < _rows; i++)
		Out.at[i][i] = 1;
	
	return Out;
}


// Matrix subtraction
extern	Matrix	subMat(Matrix _A, Matrix _B) {
	
	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'subMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_A.rows, _B.cols);

	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j] = _A.at[i][j] - _B.at[i][j];

	return Out;
}


// Multiply  matrix A and matrix B  OUT=AB
extern	Matrix	multMat(Matrix _A, Matrix _B) {
	
	if (_A.cols != _B.rows) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'multMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_A.rows, _B.cols);

	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _B.cols; j++) {
			double sum = 0;
			for (int k = 0; k < _B.rows; k++) {
				sum += _A.at[i][k] * _B.at[k][j];
			}
			Out.at[i][j] = sum;
		}
	}	
	return Out;
}


// Multiply  matrix A with a scalar k
extern	Matrix	smultMat(Matrix _A, double _k) {
	Matrix Out = createMat(_A.rows, _A.cols);
	
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			Out.at[i][j] = _k * _A.at[i][j];

	return Out;
}


// Create Transpose matrix
extern	Matrix	transpose(Matrix _A) {
	Matrix Out = createMat(_A.cols, _A.rows);
	
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			Out.at[j][i] = _A.at[i][j];

	return Out;
}


// Copy matrix
extern	Matrix	copyMat(Matrix _A) {
	Matrix Out = createMat(_A.rows, _A.cols);

	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			Out.at[i][j] = _A.at[i][j];

	return Out;

}



void gaussElim(Matrix A, Matrix b, Matrix& U, Matrix& d) {

	if (A.rows != A.cols) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: SQUARE MATRIX IS NOT UESD");
		printf("\n****************************************************\n");
	}

	int m = A.rows;
	int n = A.cols;
	double mult = 0;

	U = copyMat(A);
	d = copyMat(b);

	for (int k = 0; k < m - 1; k++) {
		for (int i = k + 1; i < m; i++) {
			mult = U.at[i][k] / U.at[k][k];
			for (int j = k; j < n; j++) {
				U.at[i][j] = U.at[i][j] - mult * U.at[k][j];
				
			}
			d.at[i][0] = d.at[i][0] - mult * d.at[k][0];
		}
	}
}


void backsub(Matrix U, Matrix& y, Matrix& x) { //Ux=y
	int n = y.rows;

	for (int i = n - 1; i >= 0; i--) {
		double sum = 0;
		for (int j = i + 1; j < n; j++) {
			sum += U.at[i][j] * x.at[j][0];
		}
		x.at[i][0] = (y.at[i][0] - sum) / U.at[i][i];
	}
}

void fwdsub(Matrix& L, Matrix& b, Matrix& y) {  //Ly=b*
	int n = b.rows;

	for (int i = 0; i < n; i++) {
		double sum = 0;
		for (int j = 0; j < i; j++) {
			sum += L.at[i][j] * y.at[j][0];
		}
		y.at[i][0] = (b.at[i][0] - sum) / L.at[i][i];
	}
}

void gaussElim_pivot(Matrix A, Matrix b, Matrix& U, Matrix& P, Matrix& L, Matrix& d) {

	if (A.rows != A.cols) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: SQUARE MATRIX IS NOT UESD");
		printf("\n****************************************************\n");
	}

	int m = A.rows;
	int n = A.cols;
	double mult = 0;
	double save_row;

	U = copyMat(A);
	d = copyMat(b);

	L = eye(m, n);
	P = eye(m, n);
	

	for (int k = 0; k < m - 1; k++) {
		int pivot = k;
		double max_value = abs(U.at[k][k]);

		for (int i = k + 1; i < m; i++) {
			if (max_value < abs(U.at[i][k])) {
				max_value = abs(U.at[i][k]);
				pivot = i;
			}
		}

		if (k != pivot) {
			for (int j = 0; j < n; j++) {
				save_row = U.at[k][j];
				U.at[k][j] = U.at[pivot][j];
				U.at[pivot][j] = save_row;
			}
			for (int j = 0; j < n; j++) {
				save_row = P.at[k][j];
				P.at[k][j] = P.at[pivot][j];
				P.at[pivot][j] = save_row;
			}
			for (int j = 0; j < k; j++) {
				save_row = L.at[k][j];
				L.at[k][j] = L.at[pivot][j];
				L.at[pivot][j] = save_row;
				
			}
			
			save_row = d.at[k][0];
			d.at[k][0] = d.at[pivot][0];
			d.at[pivot][0] = save_row;

		}

		

		for (int i = k + 1; i < m; i++) {
			mult = U.at[i][k] / U.at[k][k];
			L.at[i][k] = mult;

			for (int j = k; j < n; j++) {
				U.at[i][j] = U.at[i][j] - mult * U.at[k][j];
			}
			d.at[i][0] = d.at[i][0] - mult * d.at[k][0];
		}
	}
}



void LUdecomp(Matrix A, Matrix& L, Matrix& U, Matrix& P) {

	if (A.rows != A.cols) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: SQUARE MATRIX IS NOT UESD");
		printf("\n****************************************************\n");
	}
						//m = n
	int m = A.rows;
	int n = A.cols;
	double mult = 0;
	double save_row;

	U = copyMat(A);

	L = eye(m, n);
	P = eye(m, n);


	for (int k = 0; k < m - 1; k++) {
		Matrix Find_Max = createMat(m, 1);
		
		for (int i = k; i < n; i++) {//-------------------------------------------------------
			double Scaled_pivot_max_value = abs(U.at[i][k]); //init
			for (int j = k; j < n; j++) {
				if (Scaled_pivot_max_value < abs(U.at[i][j])) {
					Scaled_pivot_max_value = abs(U.at[i][j]);
				}
			}
			double find_max = abs(U.at[i][k]) / Scaled_pivot_max_value;
			Find_Max.at[i][0] = find_max;
		}
		int pivot = k;
		double max_value = Find_Max.at[k][0];  //init

		for (int i = n - 1; i >= 0; i--) {
			
			if (max_value <= Find_Max.at[i][0]) {
				max_value = Find_Max.at[i][0];
				pivot = i;
			}
		}
		//-----------------------------------------------------------------------------------------


		if (pivot != k) {
			for (int j = 0; j < n; j++) {
				save_row = U.at[k][j];
				U.at[k][j] = U.at[pivot][j];
				U.at[pivot][j] = save_row;
			}
			for (int j = 0; j < n; j++) {
				save_row = P.at[k][j];
				P.at[k][j] = P.at[pivot][j];
				P.at[pivot][j] = save_row;
			}
			for (int j = 0; j < k; j++) {
				save_row = L.at[k][j];
				L.at[k][j] = L.at[pivot][j];
				L.at[pivot][j] = save_row;
			}
		}


		for (int i = k + 1; i < m; i++) {
			mult = U.at[i][k] / U.at[k][k];
			L.at[i][k] = mult;

			for (int j = k; j < n; j++) {
				U.at[i][j] = U.at[i][j] - mult * U.at[k][j];
			}
		}

		//printf("At k = %d\n", k);
		//printMat(P,"[P]");
		//printMat(L, "[L_orig]");
		//printMat(U, "[U]");
		freeMat(Find_Max);
	}
}

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


double invMat(Matrix A, Matrix& Ainv) {

	int m = A.rows;
	int n = A.cols;

	if (m != n) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: SQUARE MATRIX IS NOT UESD");
		printf("\n****************************************************\n");

		return 0;
	}
	for (int i = 0; i < n; i++) {
		if (A.at[i][i] == 0) {
			printf("\n****************************************************");
			printf("\n  ERROR!!: There is zero in diagonal terms ");
			printf("\n****************************************************\n");

			return 0;
		}
	}

	Matrix L = zeros(n, n);
	Matrix U = createMat(n, n);
	Matrix P = eye(n, n);
	Matrix bb = eye(n, n);
	Matrix b = createMat(n, 1);
	Matrix x = createMat(n, 1);
	Ainv = copyMat(A);


	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			b.at[j][0] = bb.at[j][i];    //make vector cols
		}

		LUdecomp(A, L, U, P);
		solveLU(L, U, P, b, x);
		
		for (int k = 0; k < n; k++) {
			Ainv.at[k][i] = x.at[k][0];
		}
	}
	freeMat(L);
	freeMat(U);
	freeMat(P);
	freeMat(bb);
	freeMat(b);
	freeMat(x);
	return 0;
}

/////////////////////////////////////////eigvalue///////////////////////////////////////////
Matrix eigval(Matrix& A) {
	
	int m = A.rows;
	int n = A.cols;

	if (m != n) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: Not a square matrix! (A.rows != A.cols)");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	int N = 100;

	Matrix lamda = createMat(n, 1);
	Matrix U = copyMat(A);

	Matrix Q = eye(m, n);
	Matrix R = copyMat(A);

	for (int i = 0; i < N; i++) {
		QRdecomp(U, Q, R);  //QR=A
		U = multMat(R, Q);  //A=RQ 
	}
	for (int i = 0; i < n; i++) {
		lamda.at[i][0] = U.at[i][i];
	}

	//free Matrix
	freeMat(U);
	freeMat(Q);
	freeMat(R);

	return lamda;
}

void QRdecomp(Matrix& A, Matrix& Q, Matrix& R) {

	int n = A.rows;
	Matrix I = eye(n, n); //Matrix A : A.rows = A.cols
	R = copyMat(A);
	Q = copyMat(I);

	Matrix c = createMat(n, 1);
	
	for (int j = 0; j < n - 1; j++) {
		for (int k = 0; k < n; k++) {
			c.at[k][0] = R.at[k][j];
		}
		for (int k = 0; k < j; k++) {
			c.at[k][0] = 0;
		}
		Matrix e = zeros(n, 1);
		e.at[j][0] = 1;
		if (c.at[j][0] < 0) {
			e.at[j][0] = -1;
		}

		Matrix v = createMat(n, 1);
		v = addMat(c, smultMat(e, MatrixNorm(c)));  //v = c+||c||*e

		Matrix deno = multMat(transpose(v), v);     //deno = v'v

		if (deno.at[0][0] == 0) {
			printf("\n*************************************************");
			printf("\n  ERROR!!: Division by zero");
			printf("\n*************************************************\n");
			
		}
		double reciprocal = 1 / deno.at[0][0];      //reciprocal = 1/v'v
		Matrix numerator = multMat(v, transpose(v));  //numerator =vv'

		Matrix fraction = smultMat(numerator, 2 * reciprocal);  //fraction = 2(vv')/(v'v)
			
		Matrix H = subMat(I, fraction);  //H = I - 2(vv')/(v'v)

		Q = multMat(Q, H); //Q(n-1) = Q(n-2)H(n-1)
		R = multMat(H, R); //R(n-1) = H(n-1)R(n-2)

		//free Matrix
		freeMat(e);
		freeMat(v);
		freeMat(deno);
		freeMat(numerator);
		freeMat(fraction);
		freeMat(H);
	}

	//free Matrix
	freeMat(I);
	freeMat(c);
	
}

double MatrixNorm(Matrix& A) {
	int m = A.rows;  //m=3
	int n = A.cols;  //n=1

	double out = 0;

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {   
			out += pow(A.at[i][j], 2);
		}
	}
	out = pow(out, 0.5);
	return out;
}

Matrix eigvec(Matrix A) {
	int m = A.rows;
	int n = A.cols;

	if (m != n) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: Not a square matrix! (A.rows != A.cols)");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}
	if (m > 3) {
		printf("\n*************************************************");
		printf("\n  ERROR!!:Matrix A must 2x2 or 3x3");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}


	Matrix b = createMat(m - 1, n - 1); //2*2

	Matrix V = eye(m, n);
	Matrix Vi = createMat(m, 1);     //3*1
	Matrix vv = createMat(m - 1, 1); //2*1

	Matrix bb = createMat(m - 1, 1); //2*1

	for (int i = 0; i < n; i++) { //i=0,1,2 (if matA = 3x3 )

		Matrix U = copyMat(A);
		Matrix lamda = eigval(U); //3*1
		Matrix I_lamda = eye(m, n);

		for (int k = 0; k < n; k++) {
			I_lamda.at[k][k] = lamda.at[i][0];
		}

		
		Matrix B = subMat(A, I_lamda);    //B=A-lamda*I


		for (int j = 0; j < n; j++) {
			Vi.at[j][0] = V.at[j][i];
		}

		if (i == 0) {
			for (int p = 1; p < n; p++) { //1,2
				for (int q = 1; q < n; q++) { //1,2
					b.at[p - 1][q - 1] = B.at[p][q]; //(1,1), (1,2), (2,1), (2,2)
				}
				bb.at[p - 1][0] = -1 * B.at[p][i];
			}
			
			
		}
		else if (i == 1) {
			for (int p = 0; p < n; p += 2) { //0,2
				for (int q = 0; q < n; q += 2) { //0,2
					b.at[p / 2][q / 2] = B.at[p][q]; //(0,0), (0,2), (2,0), (2,2)
				}
				bb.at[p/2][0] = -1 * B.at[p][i];
			}
			
		}
		else if (i == 2) {
			for (int p = 0; p < n - 1; p++) { //0,1
				for (int q = 0; q < n - 1; q++) { //0,1
					b.at[p][q] = B.at[p][q]; //(0,0), (0,1), (1,0), (1,1)
				}
				bb.at[p][0] = -1 * B.at[p][i];
			}
			
		}

		Matrix binv = createMat(m - 1, n - 1);
	    invMat(b, binv);
		vv = multMat(binv, bb);
	   

		if (i == 0) {
			for (int h = 1; h < n; h++) { //1,2
				Vi.at[h][0] = vv.at[h - 1][0];
			}
		}
		else if (i == 1) {
			for (int h = 0; h < n; h += 2) { //0,2
				Vi.at[h][0] = vv.at[h / 2][0];
			}
		}
		else if (i == 2) {
			for (int h = 0; h < n - 1; h++) { //0,1
				Vi.at[h][0] = vv.at[h][0];
			}
		}

		if (MatrixNorm(Vi) == 0) {
			printf("\n*************************************************");
			printf("\n  ERROR!!: Division by zero");
			printf("\n*************************************************\n");

		}
		double reciprocal_norm = 1 / MatrixNorm(Vi);
		Vi = smultMat(Vi, reciprocal_norm);

		for (int h = 0; h < n; h++) { //0,1,2
			V.at[h][i] = Vi.at[h][0];
		}

		freeMat(binv);
		freeMat(U);
		freeMat(lamda);
		freeMat(I_lamda);
		freeMat(B);
	} 
	
	freeMat(b);
	freeMat(bb);
	freeMat(Vi);
	freeMat(vv);

	return V;

}

void eig(Matrix A, Matrix& V, Matrix& D) {
	
	int m = A.rows;
	int n = A.cols;

	Matrix D_eig = copyMat(D);

	Matrix lamda = eigval(A);
	Matrix sort_lamda = copyMat(lamda);

	double temp_lamda = 0.0;

	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			if (sort_lamda.at[i][0] < sort_lamda.at[j][0]) {
				temp_lamda = sort_lamda.at[i][0];
				sort_lamda.at[i][0] = sort_lamda.at[j][0];
				sort_lamda.at[j][0] = temp_lamda;
			}
		}
	}
	for (int i = 0; i < n; i++) {
		D_eig.at[i][i] = sort_lamda.at[i][0];
	}

	D = copyMat(D_eig);
	V = eigvec(A);

	freeMat(D_eig);
	freeMat(lamda);
	freeMat(sort_lamda);
}

