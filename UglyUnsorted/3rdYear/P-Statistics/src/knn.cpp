#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <string>
#include <sstream>
#include <algorithm>
#include <map>

#include "metric.h"



int main(int argc, char* argv[]){
	if (argc <= 2)
		return 1;
	
	int K;
	std::sscanf(argv[1], "%d", &K);

	// Holy crap 2 nested vectors ... 
	std::vector<std::pair<std::string, std::vector<double>>> training_samples;

	// importing training data ...
	for (int t=2; t<argc-1; t++){
		std::ifstream f;
		f.open(argv[t]);
				
		std::string temp_line;
		while (std::getline(f, temp_line)){

			training_samples.push_back(std::make_pair(std::string(argv[t]), std::vector<double>(0)));
			std::stringstream line_stream(temp_line);

			double a;
			while (line_stream >> a)
				training_samples.back().second.push_back(a);
		}
		f.close();
	}

	// import testing data ...
	std::ifstream f;
	f.open(argv[argc-1]);

	std::string sample_line;
	while (std::getline(f, sample_line)){
		std::stringstream sample_stream(sample_line);
		
		std::vector<double> sample;

		double a;
		while (sample_stream >> a)
			sample.push_back(a);
		
		std::vector<std::pair<std::string, double>> mapping_vec;
		for (size_t i=0; i<training_samples.size(); i++){
//			std::tuple<double, double, double> hpp = P::metric(training_samples[i].second, sample);

			mapping_vec.push_back(std::make_pair(training_samples[i].first, 
						std::get<0>(P::metric(training_samples[i].second, sample))));
		}
		
		// partial sorting so that K first elements are the smallest
		if (K >= mapping_vec.size()){
			std::cerr << "K larger than possible size, K=" << K << ", |z|=" << mapping_vec.size() << std::endl;
			return 1;
		}

		std::nth_element(mapping_vec.begin(), mapping_vec.begin() + K, mapping_vec.end(), 
				[](const std::pair<std::string, double>& a, 
						const std::pair<std::string, double>& b) -> bool {
					return a.second > b.second;
				});

		// now find the major group
		std::map<std::string, int> counting_map;
		for (size_t i=0; i<K; i++)
			if (counting_map.find(mapping_vec[i].first) != counting_map.end())
				counting_map[mapping_vec[i].first]++;
			else
				counting_map[mapping_vec[i].first] = 1;

//		for (auto it = mapping_vec.begin(); it != mapping_vec.end(); it++)
//			std::cout << (*it).first << " " << (*it).second << std::endl;
//		std::cout << std::endl;

		std::string group = (std::max_element(counting_map.begin(), counting_map.end(), 
				[](const std::pair<std::string, int>& a, 
						const std::pair<std::string, int>& b) -> bool {
					return a.second < b.second;
				}))->first;

		std::cout << group << std::endl;
	}

	return 0;
}
