/*
 * Viterbi.cpp 
 *
 * The calling syntax is:
 *   [vpath] = Viterbi(tau, J, M, dPara, pPara, piPara, pdfPara)
 *
 * This is a MEX file for MATLAB.
 *
*/


/* --- Headers and global variable declarations --- */

#include "mex.h"

#include "ViterbiImpl.h"
#include "InitData.h"
#include "dCalc.h"
#include "error.h"
#include "consts.h"
#include <time.h>
#include <math.h>
#include <fstream>
#include <iostream>


    double** StateIn;
	double** F;
	double** L; 
	double** G; 
	double*** H;
	double** L1;  
	double* N;  
	double** Norm;  
	double** d; 
	double** D;  
	double* mean_d;
	double** p;
	double* pi;  																																																																																						 																																																																																						
	double** eta;  																																																																																	 
	double** xi;  																																																																																	 
	double** alpha;
	int** maxI;
	int** maxU;
	double** pdf;
	//int* hiddenStates;

	int J, Y, tau, M;
	int Censoring, Output;
	bool LeftCensoring, RightCensoring;
    int* hiddenStates;


/* --- The gateway function --- */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    
    // - Input variable declarations - 
    int tauPara;
    int JPara;
    int MPara;
    double* dPara;
    double* pPara;
    double* piPara;
    double* pdfPara;  
    
    //- Output variable declarations - 
    int* hiddenStatesPara;
    
    
    // - Point to the input arguments - 

    // get scalar inputs
    tauPara = mxGetScalar(prhs[0]);
    JPara = mxGetScalar(prhs[1]);
    MPara = mxGetScalar(prhs[2]);
    
    // get arrays
    dPara = mxGetDoubles(prhs[3]);
    pPara = mxGetDoubles(prhs[4]);
    piPara = mxGetDoubles(prhs[5]);
    pdfPara = mxGetDoubles(prhs[6]);
                     
 
    // - Point to the output arguments - 
    
    // create the output matrices (adjust sizes) 
    plhs[0] = mxCreateNumericMatrix(1, tauPara, mxINT32_CLASS, mxREAL);
     
    // Point from output variables to output arguments
    hiddenStatesPara = (int*) mxGetData(plhs[0]);

    
    // - Call computational routine - 
    
    ViterbiImpl(tauPara, JPara, MPara, dPara, pPara, piPara, pdfPara, hiddenStatesPara);
    
}