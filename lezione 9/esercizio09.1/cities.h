#ifndef _cities_h_
#define _cities_h_

#include <vector>
#include <cmath>
#include "random.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

class cities{

public:
	cities(unsigned int);
	~cities(){};
	void randomOnACircumference(double);
	void randomInASquare(double);
	double distanceTraveled(vector<int>);
	void manualInitializer(vector<vector<double>>);
	double generateRandom();
	vector<vector<double>> getCoordinates();

private:
	vector<double> _cdn;
	vector<vector<double>> _coordinates;
	Random _gen;
	unsigned int _ncit; 

};

#endif
