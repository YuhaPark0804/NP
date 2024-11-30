/*-------------------------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : Yuha Park
Created          : 27-10-2024
Modified         : 27-10-2024
Language/ver     : C++ in MSVS2019

Description      : Assignment_LUdecomp
-------------------------------------------------------------------------------*/

#define ASGN		6		// enter your assignment number
#define EVAL		0		// [воик DO NOT EDIT !!!]


#include "../../include/myMatrix_22000287.h"
#include "../../include/myNP_22000287.h"


int main(int argc, char* argv[])
{

	/*	 [воик DO NOT EDIT !!!]   Resources file path setting for evaluation	*/
// #if _WIN64 | _WIN32
// 	std::string path = "C:/NP_Data/Assignment" + std::to_string(ASGN) + "/";

// #elif __APPLE__
// 	std::string path = "~/NP_Data/Assignment" + std::to_string(ASGN) + "/";
// #endif

	std::string path = "../../NP_Data/Assignment" + std::to_string(ASGN) + "/";
#if EVAL
	path += "eval/";
#endif

	/*==========================================================================*/
	/*					Variables declaration & initialization					*/
	/*--------------------------------------------------------------------------*/
	/*   - You can change the variable names									*/
	/*   - However, you must use the specified txt file name					*/
	/*==========================================================================*/

	// Option 1:  Read from datafile
	Matrix matA_Q1 = txt2Mat(path, "prob1_matA");
	Matrix vecb_Q1 = txt2Mat(path, "prob1_vecb");
	Matrix matA_Q2 = txt2Mat(path, "prob2_matK");
	Matrix vecb_Q2 = txt2Mat(path, "prob2_vecf");
	
	

	// Option 2:  Create a zero matrix with specific size
	Matrix matU_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);
	Matrix matd_Q1 = createMat(vecb_Q1.rows, vecb_Q1.cols);
	Matrix matx_Q1 = createMat(vecb_Q1.rows, vecb_Q1.cols);
	Matrix matL_Q1 = zeros(matA_Q1.rows, matA_Q1.cols); 
	Matrix matP_Q1 = eye(matA_Q1.rows, matA_Q1.cols);
	Matrix matAinv_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);

	Matrix matU_Q2 = createMat(matA_Q2.rows, matA_Q2.cols);
	Matrix matd_Q2 = createMat(vecb_Q2.rows, vecb_Q2.cols);
	Matrix matx_Q2 = createMat(vecb_Q2.rows, vecb_Q2.cols);
	Matrix matL_Q2 = zeros(matA_Q2.rows, matA_Q2.cols);
	Matrix matP_Q2 = eye(matA_Q2.rows, matA_Q2.cols);
	Matrix matAinv_Q2 = createMat(matA_Q2.rows, matA_Q2.cols);


	Matrix test_matA = txt2Mat(path, "test_matA"); //test
	Matrix test_matU = createMat(test_matA.rows, test_matA.cols); //test
	Matrix test_matL = zeros(test_matA.rows, test_matA.cols); //test
	Matrix test_matP = eye(test_matA.rows, test_matA.cols); //test



	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/

	/// Q1
	printMat(matA_Q1, "matA_Q1");
	printMat(vecb_Q1, "vecb_Q1");

	LUdecomp(matA_Q1, matL_Q1, matU_Q1, matP_Q1);
	solveLU(matL_Q1, matU_Q1, matP_Q1, vecb_Q1, matx_Q1);

	//invMat(matA_Q1, matAinv_Q1);            //Find inverse mat A
	//matx_Q1 = multMat(matAinv_Q1, vecb_Q1); //x= (A^-1) * b

	printMat(matU_Q1, "matU_Q1");
	printMat(matx_Q1, "matx_Q1");

	

	/// Q2 
	
	printMat(matA_Q2, "matA_Q2");
	printMat(vecb_Q2, "vecb_Q2");

	LUdecomp(matA_Q2, matL_Q2, matU_Q2, matP_Q2);
	solveLU(matL_Q2, matU_Q2, matP_Q2, vecb_Q2, matx_Q2);
	//invMat(matA_Q2, matAinv_Q2);              //Find inverse mat A
	//matx_Q2 = multMat(matAinv_Q2, vecb_Q2);   //x= (A^-1) * b

	printMat(matU_Q2, "matU_Q2");
	printMat(matx_Q2, "matx_Q2");
	

	//test
	 
	printMat(test_matA, "test_matA");

	LUdecomp(test_matA, test_matL, test_matU, test_matP);
	
	printMat(test_matP, "test_matP");
	printMat(test_matL, "test_matL");
	printMat(test_matU, "test_matU");

	




	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/
	freeMat(matA_Q1);		freeMat(vecb_Q1);
	freeMat(matA_Q2);		freeMat(vecb_Q2);

	freeMat(matU_Q1);		freeMat(matd_Q1);		freeMat(matx_Q1);
	freeMat(matU_Q2);		freeMat(matd_Q2);		freeMat(matx_Q2);

	freeMat(matL_Q1);		freeMat(matL_Q2);
	freeMat(matP_Q1);		freeMat(matP_Q2);

	freeMat(matAinv_Q1);
	freeMat(matAinv_Q2);


	//test
	freeMat(test_matA);
	freeMat(test_matL);
	freeMat(test_matU);
	freeMat(test_matP);

	// free other  created  matrices


	system("pause");
	return 0;
}