#ifndef _BoltzmannWeight_h_
#define _BoltzmannWeight_h_

#include "probDistr.h"
#include <cmath>

class BoltzmannWeight: public probDistr {
	public:
		void setParam(double T){_T = T;}
		double eval(double Energy) const{return exp(-Energy/_T);}
		double get_range() const{return 0.02;}
	private:
		double _T;
};

#endif
