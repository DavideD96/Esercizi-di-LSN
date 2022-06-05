#include "Metropolis.h"

using namespace std;

Metropolis::Metropolis(probDistr* funz, double startx, double starty, double startz){

	_prob_distr = funz;

	_acc = 0;
	_tent = 0;	
	
	_start[0] = startx;
	_start[1] = starty;
	_start[2] = startz;
	_position[0] = startx;
	_position[1] = starty;
	_position[2] = startz;

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

void Metropolis::move(int type){

	double range = _prob_distr->get_range();

	double x;
	double y;
	double z;

	_tent++;	

	if(type == 0){ // distribuzione uniforme
	
		//double range = 1.;
		x = _gen.Rannyu(_position[0]-range, _position[0]+range);
		y = _gen.Rannyu(_position[1]-range, _position[1]+range);
		z = _gen.Rannyu(_position[2]-range, _position[2]+range);

	} else if (type == 1){ //distribuzione gaussiana
		
		x = _gen.Gauss(_position[0], range);
		y = _gen.Gauss(_position[1], range);
		z = _gen.Gauss(_position[2], range);
	}

	double Nm = _prob_distr->eval(x,y,z);
	double Dm = _prob_distr->eval(_position[0], _position[1], _position[2]);

	double ar_probability;

	if (1.> Nm/Dm){
		ar_probability = Nm/Dm;
	}
	else{
		ar_probability = 1.;
	}

	double ar = _gen.Rannyu();
	
	if (ar < ar_probability){
		_position[0] = x;
		_position[1] = y;
		_position[2] = z;

	_acc++;
	}

}

double* Metropolis::get_position(){
	return _position;
}

void Metropolis::reset_acc_tent(){
	_acc = 0;
	_tent = 0;
}

double Metropolis::get_acc_rate(){
	
	return (double)_acc/(double)_tent;
}
