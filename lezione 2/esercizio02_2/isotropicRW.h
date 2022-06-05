/*Classe che permette di simulare un random walk con passo di lunghezza 
unitaria, orientato con probabilità uniforme su tutto l'angolo solido*/
//possibile sviluppo: fare un'unica classe RW con journey come membro pri-
//vato e clearJourney come funzione e scrivere cubicRW e isotropicRW come
//figlie.

#ifndef _isotropicRW_h_
#define _isotropicRW_h_

#include "random.h"
#include <fstream>
#include <iostream>
#include <cmath>

class isotropicRW{

public:
	isotropicRW(const double*, unsigned int); //punto di inizio, numero di passi
	~isotropicRW();
	void walk();
	void clearJourney();
	double** getJourney();
private:
	double** _journey; //riga: coordinata (0->x, 1->y, 2->z), colonna: passo
	double* _startingPosition;
	unsigned int _nStep;

	Random _gen; //se vuoi puoi anche metterlo puntatore, ricordati di mettere un new random() nel .cpp
};

#endif
