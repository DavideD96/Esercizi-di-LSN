#ifndef _Metropolis_h_
#define _Metropolis_h_

#include "probDistr.h"
#include "random.h"
#include <fstream>
#include <iostream>

class Metropolis{

public:
	Metropolis(probDistr*, double, double, double);
	~Metropolis();
	void move(int);
	double* get_position();
	void reset_acc_tent();
	double get_acc_rate();

private:
	probDistr* _prob_distr;
	Random _gen;
	double _start[3];
	double* _position = new double[3];
	unsigned int _acc, _tent;
};

#endif
