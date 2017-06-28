#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[]){
	if (argc != 4)
		return 0;

	std::ifstream a;
	std::ifstream b;
	std::ofstream c;

	std::cout << argv[1] << " " << argv[2] << " " << argv[3] << std::endl;

	a.open(argv[1]);
	b.open(argv[2]);
	c.open(argv[3]);

	for (int i=0; i<160; i++){
		for (int j=0; j<160; j++){
			int na, nb;
			a >> na; 
			b >> nb;
			c << (na + nb)/2 << " ";
		}
		c << std::endl;
	}

	a.close();
	b.close();
	c.close();

	return 0;
}
