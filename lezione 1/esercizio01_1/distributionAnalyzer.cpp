#include <cstdio>
#include "distributionAnalyzer.h"
#include <iostream>
#include <cmath>

using namespace std;

distributionAnalyzer::~distributionAnalyzer(){
	for (int i=0; i<_N; i++){
		delete[] _counts_per_bin[i];
	}
	delete[] _counts_per_bin;
}

distributionAnalyzer::distributionAnalyzer(double* data, unsigned int M){
	_randNumbers = data;
	_M = M;
	_extSx = 0;
	_extDx = 0;
	_period = M;
	_N = 1;
	_K = 1;
	_counts_per_bin = NULL;
}

void distributionAnalyzer::setExtrema(double sx, double dx){
	if(sx>=dx){
		cerr<<"Errore: estremi non validi."<<endl;
	}
	_extSx = sx;
	_extDx = dx;
}

void distributionAnalyzer::setPeriod(unsigned int period){
	_period = period;
	_N = _M/period;
}

void distributionAnalyzer::setNumBlock(unsigned int nBlocks){
	_N = nBlocks;
	_period = _M/nBlocks;
}
	
void distributionAnalyzer::setIntervals(unsigned int K){
	_K = K;
}

int** distributionAnalyzer::getCountsInInt(){
	return _counts_per_bin;
}

double* distributionAnalyzer::chi2(){

	if (_counts_per_bin == NULL){
		cerr<<"Dati per il calcolo del chi quadro assenti. "<<endl;
		return NULL;
	}
	double* chi2 = new double[_N];

	for(int i=0; i<_N; i++){
		chi2[i] = 0;
	}

	for(int j=0; j<_N; j++){
		for(int i=0; i<_K; i++){
			chi2[j] = chi2[j] + pow(_counts_per_bin[j][i]-static_cast<int>(_period/_K),2)*(static_cast<double>(_K)/_period); //per ora divisione fra interi	
		}
		//cout<<j<<endl;
		//cout<<chi2[j]<<endl;
		//cout<<endl;
	}
	return chi2;
}

void distributionAnalyzer::dataInIntervals(){

	double domain_selected[2]={_extSx,_extDx};
	double interval_lenght = (domain_selected[1]-domain_selected[0])/_K;
	int range_int[2] = {1, static_cast<int>(_K)}; //intervalli numerati da 1 a _K; man mano che vado avanti con la ricerca riduco questo range.
	int INT_NUMBER = 0; //numero dell'intervallo
	unsigned int k_;
	
	_counts_per_bin = new int*[_N];

	for (int i=0; i<_N; i++){
		_counts_per_bin[i] = new int[_K]();
	}

	for (int j=0; j<_N; j++){

		for (int i=j*_period; i < (j+1)*_period; i++){
			domain_selected[0] = _extSx;
			domain_selected[1] = _extDx;
			range_int[0] = 1;
			range_int[1] = _K;
			k_ = _K;
			INT_NUMBER = 0;
			while (INT_NUMBER == 0){
				if (k_ % 2 == 0){									//PARI
					if (_randNumbers[i] < (domain_selected[1]+domain_selected[0])/2.){		//se vale questa...
						domain_selected[1] = (domain_selected[1]+domain_selected[0])/2.;	//allora devo cercare nella metà di sinistra degli intervalli!
						range_int[1] = range_int[1] - k_/2;
					} else if (_randNumbers[i] == domain_selected[1]/2.+domain_selected[0]){
						INT_NUMBER = range_int[1] - k_/2 + 1; 					//se il numero è sul bordo appartiene all'intervallo a destra.
					} else {
						domain_selected[0] = (domain_selected[1]+domain_selected[0])/2.;
						range_int[0] = range_int[1] - k_/2 + 1;
					}
					k_ = k_/2;
				} else {										//DISPARI
					if (_randNumbers[i] < (domain_selected[1]+domain_selected[0])/2.){
						domain_selected[1] = (domain_selected[1]+domain_selected[0])/2.-interval_lenght/2.;
						range_int[1] = range_int[1] - k_/2 - 1;				//ricorda che k_/2 se k_ è dispari...	
						if (domain_selected[1] <= _randNumbers[i]){
							INT_NUMBER = range_int[1] + 1;
						}
					} else {
						domain_selected[0] = (domain_selected[1]+domain_selected[0])/2.+interval_lenght/2.;
						range_int[0] = range_int[1] - k_/2 + 1;	
						if (_randNumbers[i] < domain_selected[0]){
							INT_NUMBER = range_int[0] - 1;
						}
					}
					k_ = k_/2;	//ricorda che k_/2 se k_ è dispari...
				}
				if (range_int[0]==range_int[1] && INT_NUMBER==0){
					INT_NUMBER = range_int[0];
				}
			}
			_counts_per_bin[j][INT_NUMBER-1]++;
			//cout<<INT_NUMBER<<endl;
			//cout<<" "<<randNumbers[i]<<endl;
			//cout<<endl;
		}
	}

/*for (int i=0; i<100; i++){
	cout<<i<<endl;
	for (int j=0; j<100; j++){
		cout<<_counts_per_bin[i][j]<<" ";
	}
	cout<<endl;
}
cout<<"stop"<<endl;*/
}
