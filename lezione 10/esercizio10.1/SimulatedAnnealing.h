#ifndef _SimulatedAnnealing_h_
#define _SimulatedAnnealing_h_

#include "cities.h"
#include <vector>
#include "InitializePopulation.h"
#include <iostream>
#include <cmath>

class SimulatedAnnealing: public cities{

public:
	SimulatedAnnealing(unsigned int, double);
	~SimulatedAnnealing(){};
	void move();
	void setTemp(double);
	double searchTypicalEnergy(int);
	vector<int> getPath();
private:
	vector<int> _path;
	double _temp;
	double _energy;
	int _ncities;
};

#endif
