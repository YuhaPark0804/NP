/*-------------------------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : Yuha Park
Created          : 23-10-2024
Modified         : 23-10-2024
Language/ver     : C++ in MSVS2019

Description      : Assignment_GaussElim
-------------------------------------------------------------------------------*/

#define ASGN		5		// enter your assignment number
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
	Matrix test_matA = txt2Mat(path, "test_matA"); //optional

	// Option 2:  Create a zero matrix with specific size
	Matrix matU_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);
	Matrix matd_Q1 = createMat(vecb_Q1.rows, vecb_Q1.cols);
	Matrix matx_Q1 = createMat(vecb_Q1.rows, vecb_Q1.cols);

	Matrix matU_Q2 = createMat(matA_Q2.rows, matA_Q2.cols);
	Matrix matd_Q2 = createMat(vecb_Q2.rows, vecb_Q2.cols);
	Matrix matx_Q2 = createMat(vecb_Q2.rows, vecb_Q2.cols);

	Matrix matU_test = createMat(test_matA.rows, test_matA.cols); //optional
	Matrix matP_test = eye(test_matA.rows, test_matA.cols); //optional
	Matrix matL_test = createMat(test_matA.rows, test_matA.cols); //optional
	//Matrix matx_test = createMat(vecb_Q1.rows, vecb_Q1.cols); //optional
	//Matrix maty_test = createMat(vecb_Q1.rows, vecb_Q1.cols); //optional

	

	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/
	
	/*/// Q1
	printMat(matA_Q1, "matA_Q1");
	printMat(vecb_Q1, "vecb_Q1");

	gaussElim(matA_Q1, vecb_Q1, matU_Q1, matd_Q1);
	backsub(matU_Q1, matd_Q1, matx_Q1);

	printMat(matU_Q1, "matU_Q1");
	printMat(matd_Q1, "matd_Q1");
	printMat(matx_Q1, "matx_Q1");

	/// Q2
	printMat(matA_Q2, "matA_Q2");
	printMat(vecb_Q2, "vecb_Q2");

	gaussElim(matA_Q2, vecb_Q2, matU_Q2, matd_Q2);
	backsub(matU_Q2, matd_Q2, matx_Q2);

	printMat(matU_Q2, "matU_Q2");
	printMat(matd_Q2, "matd_Q2");
	printMat(matx_Q2, "matx_Q2");*/

	//optional
	printMat(test_matA, "test_matA");

	gaussElim_pivot(test_matA, matU_test, matP_test, matL_test);
	//LUdecomp(test_matA, matL_test, matU_test);
	//LUdecomp(test_matA, matL_test, matU_test, matP_test);
	
	//forwardsub(matL_test, matd_Q1, maty_test);
	//backsub(matU_Q1, matd_Q1, matx_test);


	printMat(matU_test, "matU_test");
	printMat(matL_test, "matL_test");
	printMat(matP_test, "matP_test");




	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/
	freeMat(matA_Q1);		freeMat(vecb_Q1); 
	freeMat(matA_Q2);		freeMat(vecb_Q2);

	freeMat(matU_Q1);		freeMat(matd_Q1);		freeMat(matx_Q1);
	freeMat(matU_Q2);		freeMat(matd_Q2);		freeMat(matx_Q2);

	//
	freeMat(test_matA);		
	freeMat(matU_test);
	freeMat(matP_test);
	freeMat(matL_test);

	// free other  created  matrices


	system("pause");
	return 0;
}