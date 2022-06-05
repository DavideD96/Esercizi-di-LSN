#ifndef _cubicRW_h_
#define _cubicRW_h_

#include "random.h"
#include <fstream>
#include <iostream>
#include <cmath>

class cubicRW{

public:
	cubicRW(const int*, unsigned int, unsigned int); //punto di inizio, dimensionalità del lattice, numero di passi
	~cubicRW();
	void walk(); //numero di cammini da ripetere 
	int** getJourney();
private:
	int** _journey; //riga: coordinata, colonna: passo
	int* _startingPosition;
	unsigned int _dim;
	unsigned int _nStep;

	Random _gen;
};

#endif
