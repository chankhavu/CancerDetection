#include <iostream>
#include <fstream>
#include <vector>
#include <string>


int main(int argc, char** argv){	
	std::vector<std::ifstream*> files(argc-2);
	int col_len = std::atoi(argv[argc-1]);

	for (int i=0; i<argc-2; i++)
		files[i] = new std::ifstream();

	for (int i=1; i<argc-1; i++)
		files[i-1]->open(argv[i]);

	std::vector<std::string> lines(argc-2);	
	bool not_empty = true;
	for (int i=0; i<files.size(); i++)
		not_empty = not_empty && (*files[i] >> lines[i]);

	while (not_empty) {
		std::string new_line;
		for (int i=0; i<lines.size(); i++){
			while (lines[i].length() < col_len)
				lines[i] = " " + lines[i];
			lines[i].resize(col_len);
			new_line = new_line + lines[i] + " ";		
		}
		std::cout << new_line << std::endl;

		for (int i=0; i<files.size(); i++)
			not_empty = not_empty && (*files[i] >> lines[i]);
	}

	for (int i=0; i<files.size(); i++)
		delete files[i];
	return 0;
}
