#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <complex>

#include <fftw3.h>
#include "EasyBMP/EasyBMP.h"
#include "otsu_threshold.h"


void f_fft(int n0, int n1, double** in, std::complex<double>** out){
	// allocate input and output
	fftw_complex* in_linear = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n0 * n1);
	fftw_complex* out_linear = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n0 * n1);

	// assign input data
	for (int i=0; i<n0; i++)
		for (int j=0; j<n1; j++){
			in_linear[i*n1 + j][0] = in[i][j];
			in_linear[i*n1 + j][1] = 0.0;
		}
	
	// perform fft
	fftw_plan p = fftw_plan_dft_2d(n0, n1, in_linear, out_linear, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);

	// assign output data
	for (int i=0; i<n0; i++)
		for (int j=0; j<n1; j++)
			out[i][j] = std::complex<double>(out_linear[i*n1 + j][0], out_linear[i*n1 + j][1]);		
	
	// cleaning
	fftw_destroy_plan(p);
	fftw_free(in_linear);
	fftw_free(out_linear);	

}
