/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : Yuha Park
Created          : 15-10-2024
Modified         : 15-10-2024
Language/ver     : C++ in MSVS2019

Description      : myMatrix.h
----------------------------------------------------------------*/

#ifndef		_MY_MATRIX_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_MATRIX_H

#include <iostream>
#include <string>
#include <fstream>

typedef struct { 
	double** at;
	int rows, cols;
}Matrix;


// Create Matrix with specified size
extern	Matrix	createMat(int _rows, int _cols);

// Free a memory allocated matrix
extern	void	freeMat(Matrix _A);

// Create a matrix from a text file
extern	Matrix	txt2Mat(std::string _filePath, std::string _fileName);

//// Print matrix
extern	void	printMat(Matrix _A, const char* _name);

// Matrix addition
extern	Matrix	addMat(Matrix _A, Matrix _B);


// initialization of Matrix elements
extern	void	initMat(Matrix _A, double _val);



//////////////////////////////////////////////////////////////////
/*							Tutorial							*/
//////////////////////////////////////////////////////////////////


// Create matrix of all zeros
extern	Matrix	zeros(int _rows, int _cols);

// Create matrix of all ones
extern	Matrix	ones(int _rows, int _cols);

// Create identity matrix
extern	Matrix	eye(int _rows, int _cols);


// Matrix subtraction
extern	Matrix	subMat(Matrix _A, Matrix _B);

// Multiply  matrix A and matrix B
extern	Matrix	multMat(Matrix _A, Matrix _B);

// Multiply  matrix A with a scalar k
extern	Matrix	smultMat(Matrix _A, double _k);

// Create Transpose matrix
extern	Matrix	transpose(Matrix _A);

// Copy matrix
extern	Matrix	copyMat(Matrix _A);


//////////////////////////////////////////////////////////////////
/*							Assginment							*/
//////////////////////////////////////////////////////////////////

void gaussElim(Matrix A, Matrix b, Matrix& U, Matrix& d);

void backsub(Matrix U, Matrix& y, Matrix& x);

void fwdsub(Matrix& L, Matrix& b, Matrix& y);

void gaussElim_pivot(Matrix A, Matrix b, Matrix& U, Matrix& P, Matrix& L, Matrix& d);

void LUdecomp(Matrix A, Matrix& L, Matrix& U, Matrix& P);

void solveLU(Matrix& L, Matrix& U, Matrix& P, Matrix& b, Matrix& x);

double invMat(Matrix A, Matrix& Ainv);

//
double MatrixNorm(Matrix& A);
Matrix eigval(Matrix& A);
void QRdecomp(Matrix& A, Matrix& Q, Matrix& R);
Matrix eigvec(Matrix A);
void eig(Matrix A, Matrix& V, Matrix& D);

#endif