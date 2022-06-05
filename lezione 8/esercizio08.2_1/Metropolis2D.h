#ifndef _Metropolis2D_h_
#define _Metropolis2D_h_

#include "probDistr.h"
#include "random.h"
#include <fstream>
#include <iostream>

class Metropolis2D{

public:
	Metropolis2D(probDistr*, double, double, int);
	~Metropolis2D();
	double* Try();
	int move(double, double);
	double* get_position();
	double* get_new_position();
	double get_Acc_rate();
	void reset_Acc_rate();
private:
	probDistr* _prob_distr;
	Random _gen;
	double _start_mu;
	double _start_sigma;
	double* _position = new double[2];
	double* _new_position = new double[2];
	int _accepted;
	int _attempted;
	int _type;
};

#endif
