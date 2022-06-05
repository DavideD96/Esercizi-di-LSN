#include "MCerror.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

MCerror::MCerror(double* resultsOfBlocks, const unsigned int n){
	if (n<=1){
	std::cerr<<"Numero di blocchi pari a 0 o 1: errore!"<<std::endl;
	}
	N_blocchi = n;
	_data = new double[N_blocchi];
	for (int i = 0; i<N_blocchi; i++){
		_data[i] = resultsOfBlocks[i];
		//std::cout<<_data[i]<<std::endl;
	}
	_errors = new double[N_blocchi];
	_means = new double[N_blocchi];
}

void MCerror::computeErrors(){ //NOTA BENE: il primo elemento dei due array è privo di significato
	double partial_A = 0.;
	double partial_A2 = 0.;
	partial_A += _data[0];
	partial_A2 += _data[0]*_data[0];
	_errors[0] = 0.; //non posso calcolare l'errore con una sola misura!
	_means[0] = partial_A;
	for(int i=1; i<N_blocchi; i++){
		partial_A += _data[i];
		partial_A2 += _data[i]*_data[i];
		_errors[i] = sqrt((1./((i+1)-1))*((1./(i+1))*partial_A2-pow(partial_A/(i+1),2)));
		_means[i] = partial_A/(i+1);
	}
}

void MCerror::setData(double* resultsOfBlocks, unsigned int n){
	if (n<=1){
	std::cerr<<"Numero di blocchi pari a 0 o 1: errore!"<<std::endl;
	}
	N_blocchi = n;
	delete[] _data;
	_data = new double[N_blocchi];
	for (int i = 0; i<N_blocchi; i++){
		_data[i] = resultsOfBlocks[i];
		//std::cout<<_data[i]<<std::endl;
	}
	_errors = new double[N_blocchi];
	_means = new double[N_blocchi];
}

double* MCerror::getErrors() const{
	return _errors;
}

double* MCerror::getMeans() const{
	return _means;
}
