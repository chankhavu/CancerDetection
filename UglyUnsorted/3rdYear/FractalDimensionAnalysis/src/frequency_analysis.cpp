#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <complex>

#include <fftw3.h>
#include "EasyBMP/EasyBMP.h"

#define IMAGE_WIDTH 160
#define IMAGE_HEIGHT 160

#define NO_BINARIZATION_MODE 0
#define BINARIZATION_MODE 1

int MODE = NO_BINARIZATION_MODE;


// assume that n0 and n1 are even numbers
double delta_beta(int n0, int n1, std::complex<double>** x){
	auto a = [&]() -> double {
		double r = 0.0;
		for (int i=0; i<n0; i++)
			for (int j=0; j<n1; j++){
				int k1 = (i < n0/2)?(i):(i - n0);
				int k2 = (j < n1/2)?(j):(j - n1);
				if (k1 != 0 && k2 != 0)
					r += std::log(std::abs(x[i][j]));
			}
		return r;
	};

	auto b = [&]() -> double {
		return n0 * n1;
	};

	auto c = [&]() -> double {
		double r = 0.0;
		for (int i=0; i<n0; i++)
			for (int j=0; j<n1; j++){
				int k1 = (i < n0/2)?(i):(i - n0);
				int k2 = (j < n1/2)?(j):(j - n1);
				if (k1 != 0 && k2 != 0)
					r += std::log(std::abs(x[i][j])) * std::log(k1*k1 + k2*k2);
			}
		return r;
	};

	auto d = [&]() -> double {
		double r = 0.0;
		for (int i=0; i<n0; i++)
			for (int j=0; j<n1; j++){
				int k1 = (i < n0/2)?(i):(i - n0);
				int k2 = (j < n1/2)?(j):(j - n1);
				if (k1 != 0 && k2 != 0)
					r += std::log(k1*k1 + k2*k2);
			}
		return r;
	};

	return a()*c() - d()*b();
}

double delta_d(int n0, int n1, std::complex<double>** x){
	auto a = [&]() -> double {
		double r = 0.0;
		for (int i=0; i<n0; i++)
			for (int j=0; j<n1; j++){
				int k1 = (i < n0/2)?(i):(i - n0);
				int k2 = (j < n1/2)?(j):(j - n1);
				if (k1 != 0 && k2 != 0)
					r += std::log(k1*k1 + k2*k2);
			}
		return r;
	};

	auto b = [&]() -> double {
		double r = 0.0;
		for (int i=0; i<n0; i++)
			for (int j=0; j<n1; j++){
				int k1 = (i < n0/2)?(i):(i - n0);
				int k2 = (j < n1/2)?(j):(j - n1);
				if (k1 != 0 && k2 != 0)
					r += std::log(std::abs(x[i][j]));
			}
		return r;
	};

	auto c = [&]() -> double {
		double r = 0.0;
		for (int i=0; i<n0; i++)
			for (int j=0; j<n1; j++){
				int k1 = (i < n0/2)?(i):(i - n0);
				int k2 = (j < n1/2)?(j):(j - n1);
				if (k1 != 0 && k2 != 0)
					r += std::log(k1*k1 + k2*k2) * std::log(k1*k1 + k2*k2);
			}
		return r;
	};

	auto d = [&]() -> double {
		double r = 0.0;
		for (int i=0; i<n0; i++)
			for (int j=0; j<n1; j++){
				int k1 = (i < n0/2)?(i):(i - n0);
				int k2 = (j < n1/2)?(j):(j - n1);
				if (k1 != 0 && k2 != 0)
					r += std::log(std::abs(x[i][j])) * std::log(k1*k1 + k2*k2);

			}
		return r;
	};

	return a()*c() - d()*b();
}

double delta(int n0, int n1, std::complex<double>** x){
	auto a = [&]() -> double {
		double r = 0.0;
		for (int i=0; i<n0; i++)
			for (int j=0; j<n1; j++){
				int k1 = (i < n0/2)?(i):(i - n0);
				int k2 = (j < n1/2)?(j):(j - n1);
				if (k1 != 0 && k2 != 0)
					r += std::log(k1*k1 + k2*k2);
			}
		return r;
	};

	auto b = [&]() -> double {
		return n0 * n1;
	};

	auto c = [&]() -> double {
		double r = 0.0;
		for (int i=0; i<n0; i++)
			for (int j=0; j<n1; j++){
				int k1 = (i < n0/2)?(i):(i - n0);
				int k2 = (j < n1/2)?(j):(j - n1);
				if (k1 != 0 && k2 != 0)
					r += std::log(k1*k1 + k2*k2) * std::log(k1*k1 + k2*k2);
			}
		return r;
	};

	auto d = [&]() -> double {
		double r = 0.0;
		for (int i=0; i<n0; i++)
			for (int j=0; j<n1; j++){
				int k1 = (i < n0/2)?(i):(i - n0);
				int k2 = (j < n1/2)?(j):(j - n1);
				if (k1 != 0 && k2 != 0)
					r += std::log(k1*k1 + k2*k2);
			}
		return r;
	};

	return a()*c() - d()*b();
}

void f_fft(int n0, int n1, double** in, std::complex<double>** out){

	// allocate input and output
	double* in_linear = (double*) malloc(sizeof(double) * n0 * n1);
	fftw_complex* out_linear = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n0 * n1);

	// assign input data
	for (int i=0; i<n0; i++)
		for (int j=0; j<n1; j++)
			in_linear[i*n1 + j] = in[i][j];

	// perform the fft 
	fftw_plan p = fftw_plan_dft_r2c_2d(n0, n1, in_linear, out_linear, FFTW_ESTIMATE);
	fftw_execute(p);

	// assign output data
	for (int i=0; i<n0; i++)
		for (int j=0; j<n1/2+1; j++){
			out[i][j] = std::complex<double>(out_linear[i*n1 + j][0], out_linear[i*n1 + j][1]);
			out[i][n1-j] = out[i][j];
		}
	
	// cleaning
	fftw_destroy_plan(p);
	free(in_linear);
	fftw_free(out_linear);	
}

double frequency_fractal_dim(int n0, int n1, int** image){
	// allocate memory
	double** im = new double* [n0];
	for (int i=0; i<n0; i++)
		im[i] = new double [n1];
	
	std::complex<double>** imfft = new std::complex<double>* [n0];
	for (int i=0; i<n0; i++)
		imfft[i] = new std::complex<double> [n1];

	// convert integer values into double
	for (int i=0; i<n0; i++)
		for (int j=0; j<n1; j++)
			im[i][j] = (double) image[i][j];

	std::cerr << "allocated" << std::endl;
	// fft
	f_fft(n0, n1, im, imfft);

	// beta = delta_beta / delta
	double beta = delta_beta(n0, n1, imfft) / delta(n0, n1, imfft);	
	double C = std::exp(delta_d(n0, n1, imfft) / delta(n0, n1, imfft));

	// slope beta is linearly related to fractal dimension D by D = E + (3-beta)/2, where
	// E is the topological dimension, which takes the value of 1 for profiles, 2 for 
	// surfaces, etc.
	double E = 2.0;
	double D = E + (3.0-beta)/2.0;

	// free memory
	for (int i=0; i<n0; i++)
		delete[] im[i];
	delete[] im;

	for (int i=0; i<n0; i++)
		delete[] imfft[i];
	delete[] imfft;

	// fractal dimension
	return D;
}


int main(int argc, char* argv[]){
	if (argc < 2)
		return 0;

#ifndef DEBUG

	std::string input_image = argv[1];
	std::ifstream img_text_file;
	img_text_file.open(input_image.c_str());
	
	// allocate memory
	int** gray_image = new int* [IMAGE_HEIGHT];
	for (int i=0; i<IMAGE_HEIGHT; i++)
		gray_image[i] = new int [IMAGE_WIDTH];

	for (int i=0; i<IMAGE_HEIGHT; i++)
		for (int j=0; j<IMAGE_WIDTH; j++){
			double a;
			img_text_file >> a;
			gray_image[i][j] = std::round(a);
		}

	img_text_file.close();

	if (argc == 3 && argv[2] == "-b")
		MODE = BINARIZATION_MODE;

	std::cerr << "processing " << argv[1] << " " << 
		((MODE==NO_BINARIZATION_MODE)?("without binarization"):("with binarization")) << std::endl;
		
	// free memory
	for (int i=0; i<IMAGE_HEIGHT; i++)
		delete[] gray_image[i];
	delete[] gray_image;

#else
	std::cerr << "start" << std::endl;
	return 0;

	BMP image;	
	image.ReadFromFile(argv[1]);

	int width = image.TellWidth();
	int height = image.TellHeight();
	int** gray_image = new int* [height];
	for (int i=0; i<height; i++)
		gray_image[i] = new int [width];

	std::cerr << "allocated" << std::endl;
	for (int i=0; i<height; i++)
		for (int j=0; j<width; j++){
			int r = (int) image(j, i)->Red;
			int g = (int) image(j, i)->Green;
			int b = (int) image(j, i)->Blue;

			gray_image[i][j] = (int) (r+g+b)/3;
		}

	std::cout << frequency_fractal_dim(height, width, gray_image);
	
	// free memory
	for (int i=0; i<height; i++)
		delete[] gray_image[i];
	delete[] gray_image;

#endif

	return 0;
}
