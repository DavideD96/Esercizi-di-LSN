#include "Metropolis2D.h"

using namespace std;

Metropolis2D::Metropolis2D(probDistr* funz, double start_mu, double start_sigma, int type){

	_type = type;

	_prob_distr = funz;

	_start_mu = start_mu;
	_start_sigma = start_sigma;

	_position[0] = start_mu;
	_position[1] = start_sigma;
	_new_position[0] = 0;
	_new_position[0] = 0;

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

Metropolis2D::~Metropolis2D(){
	delete[] _prob_distr;
}

double* Metropolis2D::Try(){

	double range = _prob_distr->get_range();

	double mu;
	double sigma;	

	if(_type == 0){ // distribuzione uniforme
	
		//double range = 1.;
		mu = _gen.Rannyu(_position[0]-range, _position[0]+range);
		sigma = _gen.Rannyu(_position[0]-range, _position[0]+range);

	} else if (_type == 1){ //distribuzione gaussiana
		
		mu = _gen.Gauss(_position[0], range);
		sigma = _gen.Gauss(_position[0], range);

	}

	_new_position[0] = mu;
	_new_position[1] = sigma;

	return _new_position;
}

int Metropolis2D::move(double E_new, double E_old){

	double Nm = _prob_distr->eval(E_new);
	double Dm = _prob_distr->eval(E_old);

	double ar_probability;

	if (1.> Nm/Dm){
		ar_probability = Nm/Dm;
	}
	else{
		ar_probability = 1.;
	}

	double ar = _gen.Rannyu();

	_attempted++;

	if (ar < ar_probability){
		_position[0] = _new_position[0];
		_position[1] = _new_position[1];
		_accepted++;

	return 1; //mossa accettata
	}else{
	
	return 0; //mossa rigettata
	}

}

double* Metropolis2D::get_position(){
	return _position;
}

double* Metropolis2D::get_new_position(){
	return _new_position;
}

double Metropolis2D::get_Acc_rate(){
	return (double)_accepted/(double)_attempted;
}

void Metropolis2D::reset_Acc_rate(){
	_accepted = 0;
	_attempted = 0;
}

