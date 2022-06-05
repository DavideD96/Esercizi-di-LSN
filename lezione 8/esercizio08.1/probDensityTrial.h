#ifndef _probDensityTrial_h_
#define _probDensityTrial_h_

#include "probDistr.h"
#include <cmath>

class probDensityTrial: public probDistr {
	public:
		//probDistr(void){_mu = 1.; _sigma = 1.;}
		void setParam(double mu, double sigma){_mu = mu; _sigma = sigma;}
		//~probDistr(){};
		double eval(double x) const{return pow(exp(-pow(x-_mu,2)/(2*_sigma*_sigma))+exp(-pow(x+_mu,2)/(2*_sigma*_sigma)),2);}
		double get_range() const{return 3.5;}
	private:
		double _mu;
		double _sigma;
		
};

#endif
