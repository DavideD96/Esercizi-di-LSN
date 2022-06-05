#ifndef _priceComputer_h_
#define _priceComputer_h_

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "random.h"

using namespace std;

struct prices{
	double call;
	double put;
};

class priceComputer{

public:
	priceComputer(double,double,double,double,double);
	~priceComputer(){};
	struct prices callPutPrice();
	double StStep(double,double);

private:
	double _initialPrice;
	double _finalInstant;
	double _strikePrice;
	double _rf_interest_rate;
	double _volatility;
	Random _gen;
	
};

#endif
