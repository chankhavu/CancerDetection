#include <iostream>
#include <fstream>
#include <vector>

std::vector<double> ordered_statistics(std::vector<double>& x);

double metric(std::vector<double>& x, std::vector<double>){
	return 0.0;
}

template <int d>
void pp(){
}

int main(int argc, char** argv){
	std::ifstream x_file;
	std::ifstream z_file;

	x_file.open(argv[1]);
	z_file.open(argv[2]);

	int f=2;
	pp<f> ();
	return 0;
}
