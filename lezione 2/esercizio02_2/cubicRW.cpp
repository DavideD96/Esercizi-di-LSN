#include "cubicRW.h"
#include <iomanip>

using namespace std;

cubicRW::cubicRW(const int start[], unsigned int dim, unsigned int npassi){
	_dim = dim;
	_nStep = npassi;
	_journey = new int*[_dim];
	_startingPosition = new int[_dim];
	
	for(int j=0; j<_dim; j++){
		_startingPosition[j] = start[j];
		_journey[j] = new int[_nStep+1]();
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

cubicRW::~cubicRW(){
	delete[] _startingPosition;
	for (int i=0; i<_dim; i++){
		delete[] _journey[i];
	}
	delete[] _journey;
}

/*void cubicRW::resetDistances(){
	_distances = {};
	std::cout<<_distances[5]<<std::endl;
}*/

void cubicRW::walk(){

	int direction;
	double randNumb;
	unsigned int coord;
	double sup = static_cast<double>(_dim);
	for(int i=0; i<_nStep; i++){
		randNumb = _gen.Rannyu(-sup, sup);
		if(randNumb == -sup){
			randNumb = 1-_dim;
		}
		coord = floor(abs(randNumb)); //seleziono la coordinata lungo la quale mi muovo.

		if(randNumb<0){
			for(int j=0; j<_dim; j++){
				_journey[j][i+1] = _journey[j][i];
			}
			_journey[coord][i+1] = _journey[coord][i]-1; 
		} else {
			for(int j=0; j<_dim; j++){
				_journey[j][i+1] = _journey[j][i];
			}
			_journey[coord][i+1] = _journey[coord][i]+1;
		}
	}
}

int** cubicRW::getJourney(){
	return _journey;
}
