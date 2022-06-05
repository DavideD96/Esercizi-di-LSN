#ifndef _MCintegral_h_
#define _MCintegral_h_
//NOTA: sostituisci il generatore con un puntatore a generatore, se lo ritieni utile

#include "Funzione.h"
#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>


class MCintegral{

public:
	MCintegral(Funzione*, double, double, unsigned int, unsigned int);
	void setNblocks(unsigned int);
	void setPeriod(unsigned int);
	~MCintegral();
	void unifSampling();
	void gaussSampling(double); //nota: questa funzione per ora si può applicare solo a integrali definiti in intervalli di R_+. La gaussiana è centrata nell'origine (questo è abbastanza facile da modificare).
	double* getResults() const;
	
private:
	unsigned int _nBlocks;
	unsigned int _period;
	unsigned int _nSampling;
	double _inf_ext;
	double _sup_ext;
	Funzione* _funz;
	double* _results;
	Random _gen;	
};

#endif
