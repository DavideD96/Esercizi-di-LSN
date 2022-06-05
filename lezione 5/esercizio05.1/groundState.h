#ifndef _groundState_h_
#define _groundState_h_

#include "probDistr.h"
#include <cmath>

class groundState: public probDistr {
	public:
		double eval(double x, double y, double z) const{return pow((1./sqrt(M_PI))*exp(-sqrt(x*x+y*y+z*z)),2);}
		double get_range() const{return 1.18;}
};

#endif
