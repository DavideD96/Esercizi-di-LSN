#include "isotropicRW.h"

using namespace std;

isotropicRW::isotropicRW(const double start[], unsigned int npassi){
	_nStep = npassi;
	_journey = new double*[3];
	_startingPosition = new double[3];
	
	for(int j=0; j<3; j++){
		_startingPosition[j] = start[j];
		_journey[j] = new double[_nStep+1]();
		_journey[j][0] = start[j];
	}
	
	//for(int i=0; i<_dim; i++){
	//	for(int j=0; j<_nStep; j++){
	//		cout<<_journey[i][j]<<" "; //era per controllare la matrice inizializzata
	//	}
	//	cout<<endl;
	//}

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

isotropicRW::~isotropicRW(){
	for (int i=0; i<3; i++){
		delete[] _journey[i];
	}
}

void isotropicRW::clearJourney(){
	for (int i=0; i<3; i++){
		for (int j=0; j<_nStep+1; j++){
		_journey[i][j] = 0.;
		}
	}
	for (int i=0; i<3; i++){
		_journey[i][0] = _startingPosition[i];
	}
}

void isotropicRW::walk(){

	double theta;
	double phi;
	for(int i = 0; i < _nStep; i++){
		theta = _gen.ZenitAngle();
		phi = _gen.Rannyu(0, 2.*M_PI);
		_journey[0][i+1] = _journey[0][i] + sin(theta)*cos(phi); 	//x = r*sin(theta)cos(phi)
		_journey[1][i+1] = _journey[1][i] + sin(theta)*sin(phi); 	//y = r*sin(theta)sin(phi)
		_journey[2][i+1] = _journey[2][i] + cos(theta);			//z = r*cos(theta)
	}
}

double** isotropicRW::getJourney(){
	return _journey;
}
