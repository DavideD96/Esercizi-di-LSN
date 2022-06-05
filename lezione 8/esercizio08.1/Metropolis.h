#ifndef _Metropolis_h_
#define _Metropolis_h_

#include "probDistr.h"
#include "random.h"
#include <fstream>
#include <iostream>

class Metropolis{

public:
	Metropolis(probDistr*, double, int);
	~Metropolis();
	void move();
	double get_position();
	double get_Acc_rate();
	void reset_Acc_rate();
private:
	probDistr* _prob_distr;
	Random _gen;
	double _start;
	double _position;
	int _accepted;
	int _attempted;
	int _type;
};

#endif
