#include <fstream>
#include <iostream>

int main(int argc, char** argv){
	std::ifstream f;
	f.open(argv[1]);

	double x, min = 1000.0, max = -1000.0;
	while (f >> x){
		min = (x < min)?(x):(min);
		max = (x > max)?(x):(max);
	}

	std::cout << "min element: " << min << std::endl;
	std::cout << "max element: " << max << std::endl;
	return 0;
}
