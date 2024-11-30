/*------------------------------------------------------------------------------------------\
@ Numerical Programming  by Young-Keun Kim - Handong Global University

Author          : Young-Keun Kim
Created         : 01-04-2019
Modified        : 01-04-2023
Language/ver	: C in MSVS2017
Course			: Numerical Programming

Description     : Assignment 7 Eigenvalue,Eigenvector
/------------------------------------------------------------------------------------------*/

#define ASGN	7		// enter your assignment number
#define EVAL	0		// set 0

// [ MODIFY Path and name of Headerfiles]
#include "../../include/myMatrix_22000287.h"
#include "../../include/myNP_22000287.h"


// [ CREATE FUNCTIONS ]
// These must be in myMatrix.h, myMatrix.cpp
// Matrix eigval(Matrix A);
// void QRdecomp(Matrix A, Matrix Q, Matrix R);
// Matrix eigvec(Matrix A);
// void eig(Matrix A, Matrix V, Matrix D);




int main(int argc, char* argv[])
{
#if _WIN64 | _WIN32
	/*	 [���� DO NOT EDIT !!!]   Resources file path setting for evaluation	*/
	std::string path = "../../NP_Data/Assignment" + std::to_string(ASGN) + "/";
#elif __APPLE__
	//std::string path = "~/NP_Data/Assignment" + std::to_string(ASGN) + "/";
	std::string path = "../../NP_Data/Assignment" + std::to_string(ASGN) + "/";
#endif
#if EVAL
	path += "eval/";
#endif

	/*==========================================================================*/
	/*					Variables declaration & initialization					*/
	/*--------------------------------------------------------------------------*/
	/*   - You can change the variable names									*/
	/*   - However, you must use the specified file name						*/
	/*		   : For each assignment, the file name will be notified			*/
	/*==========================================================================*/
	Matrix matA_Q1 = txt2Mat(path, "prob1_matA");
	Matrix eigVals_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);
	Matrix eigVecs_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);
	Matrix matD_Q1 = eye(matA_Q1.rows, matA_Q1.cols);
	Matrix matV_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);


	Matrix matA_Q2 = txt2Mat(path, "prob2_matA");
	Matrix eigVals_Q2 = createMat(matA_Q2.rows, matA_Q2.cols);
	Matrix eigVecs_Q2 = createMat(matA_Q2.rows, matA_Q2.cols);
	Matrix matD_Q2 = eye(matA_Q2.rows, matA_Q2.cols);
	Matrix matV_Q2 = createMat(matA_Q2.rows, matA_Q2.cols);


	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/

	// [Eigenval/vectors for Q1]
	
	//eigVals_Q1 = eigval(matA_Q1);
	//eigVecs_Q1 = eigvec(matA_Q1);
	eig(matA_Q1, matV_Q1, matD_Q1);


	// [Eigenval/vectors for Q2]
	
	//eigVals_Q2 = eigval(matA_Q2);
	//eigVecs_Q2 = eigvec(matA_Q2);
	eig(matA_Q2, matV_Q2, matD_Q2);


	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/

	// [Eigenval/vectors for Q1]
	printf("----------------------------------------------------------------------------------------------\n");
	printf("			       Eigenvalue & Eigenvector Results					\n");
	printf("----------------------------------------------------------------------------------------------\n");
	printf("\n[Eigen value]\n\n");
	//printMat(eigVals_Q1,"eigVals_Q1");
	printMat(matD_Q1, "matD_Q1");

	printf("\n[Eigen vector]\n\n");
	//printMat(eigVecs_Q1,"eigVecs_Q1");
	printMat(matV_Q1, "matV_Q1");
	printf("\n\n");


	// [Eigenval/vectors for Q2]
	printf("\n[Eigen value]\n\n");
	//printMat(eigVals_Q2, "eigVals_Q2");
	printMat(matD_Q2, "matD_Q2");

	printf("\n[Eigen vector]\n\n");
	//printMat(eigVecs_Q2, "eigVecs_Q2");
	printMat(matV_Q2, "matV_Q2");
	printf("\n\n");




	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/

	// [Deallocate for Q1]
	freeMat(matA_Q1);
	freeMat(eigVals_Q1);	freeMat(eigVecs_Q1);
	freeMat(matD_Q1);       freeMat(matV_Q1);


	// [Deallocate for Q2]
	freeMat(matA_Q2);
	freeMat(eigVals_Q2);	freeMat(eigVecs_Q2);
	freeMat(matD_Q2);       freeMat(matV_Q2);


	system("pause");
	return 0;
}