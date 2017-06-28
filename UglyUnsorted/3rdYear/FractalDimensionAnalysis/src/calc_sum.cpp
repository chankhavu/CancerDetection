#include <fstream>
#include <iostream>

int main(int argc, char** argv){
	std::ifstream f;
	f.open(argv[1]);
	double s = 0.0, x;
	int c = 0;
	while (f >> x){
		c++;
		s += x;
	}

	f.close();
	std::cout << s / (1.0 * c);
}
