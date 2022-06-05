#ifndef _MCerror_h_
#define _MCerror_h_

class MCerror {

public:
	MCerror(double*, unsigned int); //costruttore: inserire il numero di blocchi per la simulazione
	~MCerror(){delete[] _data;};
	void setData(double*, unsigned int);
	void computeErrors();
	double* getErrors() const;
	double* getMeans() const;
private:
	double* _data;
	unsigned int N_blocchi;
	double* _errors;
	double* _means;
};

#endif
