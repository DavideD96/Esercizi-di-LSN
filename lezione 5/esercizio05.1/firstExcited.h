#ifndef _firstExcited_h_
#define _firstExcited_h_

#include "probDistr.h"
#include <cmath>

class firstExcited: public probDistr {
	public:
		double eval(double x, double y, double z) const{return pow((1./(8*sqrt(2./M_PI)))*sqrt(x*x+y*y+z*z)*exp(-sqrt(x*x+y*y+z*z)/2)*z/sqrt(x*x+y*y+z*z),2);}
		double get_range() const{return 2.8;}
};

#endif
