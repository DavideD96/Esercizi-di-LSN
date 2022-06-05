#ifndef _geneticAlgorithm_h
#define _geneticAlgorithm_h

#include "cities.h"
#include <vector>
#include "InitializePopulation.h"
#include <iostream>
#include <cmath>

using namespace std;

class geneticAlgorithm: public cities{

public:
	geneticAlgorithm(unsigned int, unsigned int);
	~geneticAlgorithm(){};
	void initializePop(int);
	void selection();
	vector<vector<int>> getPopulation();
	double meanBestHalf();
	double minLenght();
	void mutation_permutation(); //simple permutation
	void mutation_inversion_mcities(); //inverts the order of a sequence of m cities
	void mutation_exchange(); //exchange two sequences of m cities
	void crossover();
	int find(vector<int>,int);
	int findMin(vector<int>);
	vector<int> bestTravel();
	void replacePath(vector<int>);
	void pointReplace(int,int,int); //function for debug: n to replace, path, step.

private:
	vector<int> _singlePath;
	vector<double> _distances;
	vector<vector<int>> _population;
	unsigned int _ncities;
	unsigned int _npopulation;
};

#endif
