#ifndef _cosine_h_
#define _cosine_h_

#include "Funzione.h"
#include <cmath>

class cosine: public Funzione {
	public:
		double eval(double x) const{return (M_PI/2.)*cos((M_PI/2.)*x);}
};

#endif
		
