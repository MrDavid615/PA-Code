#include <stdio.h>
#include <mex.h>
#include <math.h>
#include <iostream>

// C”Ô—‘∞¥––±È¿˙

void FFBP_fun_C(int pixel_Num, double Signal[][128] , int POINT[][512*512], double COS[][512*512], double (*ROI)[2])
{
	for (int t = 0; t < pixel_Num; t++)
	{
		for (int i = 0; i < 128; i++)
		{
			ROI[t][0] -= Signal[(POINT[i][t])][i] * COS[i][t];
		}
        ROI[t][1] = -ROI[t][0];
	}
}


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	int pixel_Num = mxGetScalar(prhs[0]);
	double (*Signal)[128] = (double (*)[128])mxGetPr(prhs[1]);
	int (*POINT)[512*512] = (int (*)[512*512])mxGetPr(prhs[2]);
	double (*COS)[512*512] = (double(*)[512*512])mxGetPr(prhs[3]);

	plhs[0] = mxCreateDoubleMatrix(2, pixel_Num, mxREAL);
	double (*ROI)[2];
	ROI = (double (*)[2])mxGetPr(plhs[0]);
	FFBP_fun_C(pixel_Num, Signal, POINT, COS, ROI);
}



