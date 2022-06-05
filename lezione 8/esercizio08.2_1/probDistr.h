#ifndef _probDistr_h_
#define _probDistr_h_

class probDistr{
	public: 
		virtual double eval(double x) const = 0;
		virtual double get_range() const = 0;
};

#endif
