#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

/**
* Custom struct to store the evaluation of a sequence.
* @param cs_score_ CS-Score of a sequence
* @param cs_percentage_ percentage of residues above threshold (CS-Score)
* @param css_score_ CSS_Score of a sequence
* @param css_percentage_ percentage of residues above threshold (CSS-Score)
* @param sequence_number_ id of sequence 
*/

struct Entry{
	double cs_score_;
	double cs_percentage_;
	double css_score_;
	double css_percentage_;
	int sequence_number_;
};
/**
* Compares two entries
* @param a one entry of type Entry
* @param b second entry of type Entry
* @return true if CS percentage of a is bigger than b otherwise false
*/
bool comparator(Entry a, Entry b){
	return a.cs_percentage_ > b.cs_percentage_;
}

/**
* Extracts value from string
* @param content string from which value should be extracted
* @return value extracted from string
*/
double getValue(std::string& content){
	std::string delimiter = " ";
	std::vector<double> str;
	auto pos = 0;
	
	while((pos = content.find(delimiter)) != std::string::npos){
		std::string s = content.substr(0, pos);
		try{
			double x = std::stod(s);
			str.push_back(x);
			content.erase(0, pos + delimiter.length());
		}
		catch (const std::invalid_argument& ia){
			content.erase(0, pos + delimiter.length());
			continue;
			
		}
	}
	
	try{
		double x = std::stod(content);
		str.push_back(x);
	}catch (const std::invalid_argument& ia){
	}
	
	
	return str[0];
}



int main(int argc, char* argv[]) {
    if(argc != 2){
        std::cerr << "Wrong number of arguments given. Should be 1 but was " + (argc - 1) << std::endl;
	exit(1);
    }
	
    auto output = argv[1];
	
	std::string folder = "../Evaluation/";
	std::vector<Entry> sorting;
	
	for(auto i = 0; i < 999; i++){
		std::string file = folder + "Sequence_" + std::to_string(i) + ".DAT" ;
		
		try{
			std::ifstream in;
			in.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
			in.open(file);
		}catch (const std::ifstream::failure& e) {
			std::cout << "Could not process sequence number " << i << std::endl;
			continue;
        	}
		std::ifstream in(file);
		std::vector<std::string> content; 
		std::string line;
		std::string delimiter = ":";
		int j = 0;
		
		while(std::getline(in, line)){
			auto pos = line.find(delimiter);
			if(pos != std::string::npos){
				line.erase(0, pos + delimiter.length());
				content.push_back(line);
			}
		}
		
		content.erase(content.begin());
		
		std::vector<double> info; 
		
		for(auto e:content){
			double s = getValue(e);
			info.push_back(s);
		}
		
		Entry e;
		e.cs_score_ = info[0];
		e.cs_percentage_ = info[1];
		e.css_score_ = info[2];
		e.css_percentage_ = info[3];
		e.sequence_number_ = i;
		
		sorting.push_back(e);	
	}
 
	std::sort(sorting.begin(), sorting.end(), comparator);
	
	std::ofstream out(output);
	
	out << "Sequence number\tCS-Score\tCS-Percentage\tCSS-Score\tCSS-Percentage\t" << std::endl;  
	
	for(auto e:sorting){
		if(e.cs_score_ < 32.15){
			continue;
		}else{
			out << e.sequence_number_ << "\t" << e.cs_score_ << "\t" << e.cs_percentage_ << "\t" << e.css_score_ << "\t" << e.css_percentage_ << "\t" << std::endl;
		}
		
	}

  return 0;
}

