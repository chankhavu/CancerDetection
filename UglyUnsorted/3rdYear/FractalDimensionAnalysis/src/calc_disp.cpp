#include <fstream>
#include <sstream>
#include <iostream>

int main(int argc, char** argv){
	std::ifstream f1, f2;
	f1.open(argv[1]);
	f2.open(argv[2]);

	std::string s1, s2;

	// f1 and f2 are designed so that the number of lines are 
	// equal. A line can be empty if the matching bmp image
	// is damaged.
	while (std::getline(f1, s1) && std::getline(f2, s2)){
		if (s1 == "" || s2 == ""){
			std::cerr << "empty line found in files " << argv[1] << ", " << argv[2] << std::endl;
			continue;
		}

		std::stringstream ss1(s1), ss2(s2);
		double n1, n2;
		ss1 >> n1;
		ss2 >> n2;

		std::cout << n2-n1 << std::endl;
	}

	return 0;
}
