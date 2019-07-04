/*
 * FB.c 
 *
 * The calling syntax is:
 *   [F,L,G,L1,N,Norm,eta,xi,error] = FB(censoring, tau, J, M, 
 *                                        dPara, pPara, piPara, pdfPara, 
 *                                          F, L, G, L1, N, Norm, eta, xi, err)
 *
 * This is a MEX file for MATLAB.
 *
*/


/* --- Headers and global variable declarations --- */

#include "mex.h"

#include "FBImpl.h"
#include "InitData.h"
#include "dCalc.h"
#include "error.h"
#include "mymat.h"
#include "cube.h"
#include "consts.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <limits.h>
#include <iostream>
#include <fstream>


/* --- The gateway function --- */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    
    // - Input variable declarations - 
    int censoringPara;
    int tauPara;
    int JPara;
    int MPara;
    double* dPara;
    double* pPara;
    double* piPara;
    double* pdfPara;  
    
    //- Output variable declarations - 
    double* FPara;
    double* LPara;
    double* GPara;
    double* L1Para;
    double* NPara;
    double* NormPara;
    double* etaPara;
    double* xiPara;
    int errorPara;
    
    
    // - Point to the input arguments - 

    // get scalar inputs
    censoringPara = mxGetScalar(prhs[0]); 
    tauPara = mxGetScalar(prhs[1]);
    JPara = mxGetScalar(prhs[2]);
    MPara = mxGetScalar(prhs[3]);
    
    errorPara = mxGetScalar(prhs[16]);
            
    // get arrays
    dPara = mxGetDoubles(prhs[4]);
    pPara = mxGetDoubles(prhs[5]);
    piPara = mxGetDoubles(prhs[6]);
    pdfPara = mxGetDoubles(prhs[7]);
                                 
    // - Point to the output arguments - 
    
    // create the output matrices (adjust sizes) 
    
    plhs[0] = mxCreateDoubleMatrix(1, JPara * tauPara, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, JPara * tauPara, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, JPara * tauPara, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1, JPara * tauPara, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1, tauPara, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(1, JPara * tauPara, mxREAL);
    plhs[6] = mxCreateDoubleMatrix(1, JPara * MPara, mxREAL);
    plhs[7] = mxCreateDoubleMatrix(1, JPara * MPara, mxREAL);
    plhs[8] = mxCreateDoubleMatrix(1, 1, mxREAL);
     
    // Point from output variables to output arguments
    FPara = mxGetDoubles(plhs[0]);
    LPara = mxGetDoubles(plhs[1]);
    GPara = mxGetDoubles(plhs[2]);
    L1Para = mxGetDoubles(plhs[3]);
    NPara = mxGetDoubles(plhs[4]);
    NormPara = mxGetDoubles(plhs[5]);
    etaPara = mxGetDoubles(plhs[6]);
    xiPara = mxGetDoubles(plhs[7]);
    errorPara = mxGetScalar(plhs[8]);
    
    // - Call computational routine - 
    
    FBImpl(censoringPara, tauPara, JPara, MPara, dPara, pPara, piPara, pdfPara, FPara, LPara, GPara, L1Para, NPara, NormPara, etaPara, xiPara, &errorPara);

}