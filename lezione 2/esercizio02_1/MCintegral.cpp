#include "MCintegral.h"

using namespace std;

MCintegral::MCintegral(Funzione* funz, double a, double b, unsigned int N, unsigned int P){
	_nBlocks = N;
	_period = P;
	_nSampling = N*P;
	_inf_ext = a;
	_sup_ext = b;
	_funz = funz;
	_results = new double[_nBlocks];

	ifstream inData;

	int seed[4];
	int p1, p2;

	inData.open("Primes.txt");

	if (inData.is_open()){
		inData >> p1 >> p2;
	} 
	else{
		std::cerr << "Unable to open Primes" <<std::endl;
	}
	inData.close(); 

	inData.open("seed.in");
	string property;

	if (inData.is_open()){
		while(!inData.eof()){
			inData >> property;
			if(property == "RANDOMSEED"){
				inData >> seed[0] >> seed[1] >> seed[2] >>seed[3];
			_gen.SetRandom(seed, p1, p2);
			}
		}
		inData.close();
	}

}

MCintegral::~MCintegral(){
	delete[] _funz;
	delete[] _results;
}

void MCintegral::setNblocks(unsigned int N){
	_nBlocks = N;
	_nSampling = _nBlocks*_period;
	delete[] _results;
	_results = new double[_nBlocks];
}

void MCintegral::setPeriod(unsigned int P){
	_period = P;
	_nSampling = _nBlocks*_period;
}

void MCintegral::unifSampling(){
	double a;
	double sumPart;
	for(int j=0; j<_nBlocks; j++){
		sumPart = 0.;
		for(int i=0; i<_period; i++){
			a = _gen.Rannyu(_inf_ext, _sup_ext);
			sumPart+=_funz->eval(a);
		}
		_results[j]=(_sup_ext-_inf_ext)*sumPart/_period;
		//std::cout<<_results[j]<<std::endl;
	}
}

void MCintegral::gaussSampling(double sigma){
	double mean = 0;
	int i;
	double singleRand;
	double sumPart;
	double denom;
	for(int j=0; j<_nBlocks; j++){
		sumPart = 0.;
		i = 0;
		while(i<_period){
			singleRand = _gen.Gauss(mean, sigma);
			if (singleRand>=_inf_ext && singleRand<_sup_ext){
				//cout<<singleRand<<endl;
				denom = (1./(erf(_sup_ext/(sigma*sqrt(2)))-erf(_inf_ext/(sigma*sqrt(2)))))*(2./sqrt(2*M_PI*pow(sigma,2)))*(exp(-0.5*pow((singleRand-mean)/sigma,2)));
				sumPart += _funz->eval(singleRand)/denom;
				i++;
			}
		}
		_results[j] = (_sup_ext-_inf_ext)*sumPart/_period;
		//cout<<_results[j]<<endl;
		//cout<<erf(_sup_ext/(sigma*sqrt(2)))<<endl;
	}
}

double* MCintegral::getResults() const{
	return _results;
}
