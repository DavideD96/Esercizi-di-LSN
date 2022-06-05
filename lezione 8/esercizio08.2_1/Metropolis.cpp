#include "Metropolis.h"

using namespace std;

Metropolis::Metropolis(probDistr* funz, double startx, int type){

	_prob_distr = funz;

	_type = type;

	_start = startx;

	_position = startx;

	_accepted = 0;
	_attempted = 0;

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

Metropolis::~Metropolis(){
	delete[] _prob_distr;
}

void Metropolis::move(){

	double range = _prob_distr->get_range();

	double x;	

	if(_type == 0){ // distribuzione uniforme
	
		//double range = 1.;
		x = _gen.Rannyu(_position-range, _position+range);

	} else if (_type == 1){ //distribuzione gaussiana
		
		x = _gen.Gauss(_position, range);

	}

	double Nm = _prob_distr->eval(x);
	double Dm = _prob_distr->eval(_position);

	double ar_probability;

	if (1.> Nm/Dm){
		ar_probability = Nm/Dm;
	}
	else{
		ar_probability = 1.;
	}

	double ar = _gen.Rannyu();
	
	if (ar < ar_probability){
		_position = x;
		_accepted++;
	}
	_attempted++;

}

double Metropolis::get_position(){
	return _position;
}

double Metropolis::get_Acc_rate(){
	return (double)_accepted/(double)_attempted;
}

void Metropolis::reset_Acc_rate(){
	_accepted = 0;
	_attempted = 0;
}

