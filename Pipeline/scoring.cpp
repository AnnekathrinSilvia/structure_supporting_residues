#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <iterator>

/**
* This script calculates the conservation scores for each position of a given input alignment. 
* The used conservation score scheme was developed by Valdar et al. (https://doi.org/10.1002/1097-0134(20010101)42:1<108::AID-PROT110>3.0.CO;2-O).
*/


/**
* Finds the largest value of a matrix
*
* @param matrix the matrix whose largest value should be found
* @return the largest value of the matrix
*/
double getMax(std::vector<std::vector<double>>& matrix){
	double max = std::numeric_limits<double>::min();
	
	for (auto it = matrix.begin(); it != matrix.end(); it++){
		for(auto it2 = it->begin(); it2 != it->end(); it2++){
			if(*it2 > max){
				max = *it2;
			}
		}
	}
	
	return max;
}

/**
* Finds the smallest value of a matrix
*
* @param matrix the matrix whose smallest value should be found
* @return the smallest value of the matrix
*/
double getMin(std::vector<std::vector<double>>& matrix){
	double min = std::numeric_limits<double>::max();
	
	for(auto it = matrix.begin(); it != matrix.end(); it++){
		for(auto it2 = it->begin(); it2 != it->end(); it2++){
			if(*it2 < min){
				min = *it2;
			}
		}
	}
	
	return min;
}

/**
* Calculates Lambda, the value to normalize the scores
* @param n number of sequences
* @param weights vector of weights of size name
* @return reciprocal of sum of products of individual weights
*/
long double calculateLambda(std::size_t n, std::vector<long double> weights){
    long double sum = 0;
    for(auto x = 0; x < n; x++){
        long double w = weights[x];
        for(auto y = x+1; y < n; y++){
            sum += w * weights[y];
        }
    }
    return 1 / sum;
}




/**
* Retriev substitution value for two amino acids
* @param pos position index of amino acids which should be substituted
* @param seq_i one of two sequences to retrieve amino acid
* @param seq_j second of two sequences to retrieve amino acid
* @param matrix substitution matrix
* @param amino_acids vector of amino acids in matrix
* @param max largest value of matrix
* @param min smallest value of matrix
* @return substitution value for two amino acids
*/
double getMatrix(std::size_t pos,std::string& seq_i, std::string& seq_j,std::vector<std::vector<double>>& matrix, std::vector<char>& amino_acids, double max, double min){
	///Retrieve amino acids
    char x = seq_i[pos];
    char y = seq_j[pos];

	///Verify that neither x or y are gaps
    if(x == '.' || y == '.' || x == '-' || y == '-'){
        return 0;
    }
	
	///Retrieve index of amino acids to be substitution
	auto index_x = std::distance(amino_acids.begin(), std::find(amino_acids.begin(), amino_acids.end(), x));
	auto index_y = std::distance(amino_acids.begin(), std::find(amino_acids.begin(), amino_acids.end(), y));
	
	///Get substitution value
    double m_ij = matrix[index_x][index_y];
	///Normalization of substitution value to ensure bounded output space of conservation values
    auto nom = m_ij - min;
    auto denom = max - min;
    return nom/denom;
}

/**
* Calculate conservation score of one position
* @param pos position for which conservation score should be calculated
* @param sequences vector of protein sequences
* @param matrix given substitution matrix
* @param amino_acids vector of amino acids in matrix
* @param weight vector of weights of length(sequences)
* @param max largest value of given substitution matrix
* @param min smallest value of given substitution matrix
* @return conservation score for pos
*/
long double calculateScore(std::size_t pos, std::vector<std::string>& sequences, std::vector<std::vector<double>>& matrix, std::vector<char>& amino_acids, std::vector<long double>& weights, double max, double min){
    long double score = 0;
    for(auto i = 0; i < sequences.size(); i++){
        auto w = weights[i];
        for(auto j = i+1; j < sequences.size(); j++){
            auto m = getMatrix(pos, sequences[i], sequences[j], matrix, amino_acids, max, min);
            score += w * weights[j] * m;
        }
    }
    return score;
}

/**
* Determines all non-gap positions of two sequences
* @param seq_i one of two protein sequence
* @param seq_j second protein sequence
* @return set of non-gap positions
*/
std::vector<int> getNoGapPositions(std::string& seq_i, std::string& seq_j){
    std::vector<int> pos;
    for(auto i = 0; i < seq_i.size(); i++){
        char x = seq_i[i];
        char y = seq_j[i];
		
		///Check if x AND y are gap
        if((x == '.' && y == '.') || (x == '-' && y == '-')){
            continue;
        }else{
	        pos.push_back(i);
        }
    }
    return pos;
}

/**
* Calculates distance of two sequences
* @param seq_i one of two sequences
* @param seq_j second of two sequences
* @param matrix given substitution matrix
* @param amino_acids amino acids in matrix
* @param max largest value of given substitution matrix
* @param min smallest value of given substitution matrix
* @return distance of seq_i and seq_j
*/
long double getDistance(std::string& seq_i, std::string& seq_j, std::vector<std::vector<double>>& matrix, std::vector<char>& amino_acids, double max,double min){
    long double sum = 0;
    auto positions = getNoGapPositions(seq_i, seq_j);
    for(auto e:positions){
        sum += getMatrix(e, seq_i, seq_j, matrix, amino_acids, max, min);
    }
    return 1 - (sum / positions.size());
}

/**
* Calculates weight of one sequence
* @param i index of considered sequence in sequence vector
* @param sequences vector of protein sequences
* @param matrix given substitution matrix
* @param amino_acids amino acids in matrix
* @param max largest value of given substitution matrix
* @param min smallest value of given substitution matrix
* @return weight of sequence
*/
long double getWeight(std::size_t i, std::vector<std::string>& sequences, std::vector<std::vector<double>>& matrix, std::vector<char>& amino_acids, double max,double min, std::vector<std::vector<long double>>& distances){
    long double sum = 0;
    for(std::size_t k = 0; k < sequences.size(); k++){
        if(k != i){
            auto w = distances[i];
            if(w[k] != std::numeric_limits<long double>::min()){
                sum += w[k];
            }else{
                auto x = getDistance(sequences[i], sequences[k], matrix, amino_acids, max, min);
                sum += x;
                w[k] = x;
                distances[i] = w;
            }
        }
    }
    return sum / (sequences.size() - 1);
}

int main(int argc, char* argv[]) {
    if(argc != 4){
        std::cerr << "Wrong number of arguments given. Should be 3 but was " + (argc - 1) << std::endl;
	exit(1);
    }

    auto input = argv[1];
    auto mat = argv[2];
	auto output = argv[3];

    std::ifstream in(input);

    std::vector<std::string> sequences;
    std::string current_sequence;
    std::string line;

	///Read sequences of alignment into vector 
    while(std::getline(in, line)){
        if(line[0] == '>'){
            sequences.push_back(current_sequence);
            current_sequence.clear();
        }else{
            current_sequence += line;
        }
    }
    sequences.erase(sequences.begin());
    sequences.push_back(current_sequence);
    std::cout << "Read sequences" << std::endl;

    ///Number of sequences in the alignment
    std::size_t n = sequences.size();

    std::ifstream m(mat);
	std::vector<std::vector<double>> matrix;
    std::vector<char> amino_acids;
    int a = 0;
    std::string str;
	
	
	///Read substutiton matrix into matrix format
    while(std::getline(m, str)){
        if(a == 0){
            a = 2;
            auto delimiter = " ";
            std::string::size_type beg = 0;

            for(auto end = 0; (end = str.find(delimiter, end)) != std::string::npos; end++){
                std::string s = str.substr(beg, end - beg);
                amino_acids.push_back(s[0]);
                beg = end + 1;
            }
        }else{
            auto delimiter = " ";
            std::string::size_type beg = 0;
            std::vector<double> s;
			int i = 0;

            for(auto end = 0; (end = str.find(delimiter, end)) != std::string::npos; end++){
				if(i == 0){
					i = 1;
					beg = end + 1;
					std::string x = str.substr(beg, end);
					continue;
				}
                std::string x = str.substr(beg, end);
                s.push_back(std::stod(x));
                beg = end + 1;
            }
            matrix.push_back(s);
        }
    }

    std::cout << "Build matrix" << std::endl;
	 
	 
    auto max = getMax(matrix);
    auto min = getMin(matrix);
	 
	std::cout << "Retrieved maximal and minimal value of the matrix" << std::endl;

    std::vector<long double> weights;
    std::vector<std::vector<long double>> distance;
    for(std::size_t i = 0; i < n; i++){
        std::vector<long double> d(n, std::numeric_limits<long double>::min());
        distance.push_back(d);
    }
    std::cout << "Build distance" << std::endl;

    for(std::size_t i = 0; i < n; i++){
        auto w = getWeight(i, sequences, matrix, amino_acids, max, min, distance);
        weights.push_back(w);
    }

    std::cout << "Calculated weights" << std::endl;

    auto lamb = calculateLambda(n, weights);

    std::cout << "Calculated Lambda" << std::endl;

    std::ofstream o;
    o.open(output);

    std::size_t len = sequences[0].size();

    for(auto i = 0; i < len; i++){
        auto score = calculateScore(i, sequences, matrix, amino_acids, weights, max, min);
        score *= lamb;
        o << score << "\n";
    }
	
  return 0;
}

