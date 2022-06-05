#include "cities.h"

cities::cities(unsigned int ncit){
	
	_ncit = ncit;
	_cdn.resize(3,0);
	
	_coordinates.resize(ncit,_cdn);

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

vector<vector<double>> cities::getCoordinates(){
	
	return _coordinates;
}

void cities::randomOnACircumference(double radius){

	double theta, x, y;

	for(int i=0; i<_ncit; i++){

		theta = _gen.Rannyu(0,2*M_PI);

		x = radius*cos(theta);
		y = radius*sin(theta);

		_coordinates[i][0] = i+1;
		_coordinates[i][1] = x;
		_coordinates[i][2] = y;

	}
	
}

void cities::randomInASquare(double lato){

	for(int i=0; i<_ncit; i++){

		_coordinates[i][0] = i+1;
		_coordinates[i][1] = _gen.Rannyu(0,lato);
		_coordinates[i][2] = _gen.Rannyu(0,lato);

	}
}

double cities::distanceTraveled(vector<int> sequence){
	int attuale, precedente;
	double distance2 = 0;	

	/*for(int i=0; i<sequence.size(); i++){
		cout << sequence[i] << endl;
	}*/

	for(int i=1; i<sequence.size(); i++){
		attuale = sequence[i]-1;
		precedente = sequence[i-1]-1;
		//cout << "x1 = "<<_coordinates[attuale][1]<< " y1 = " <<_coordinates[attuale][2]<< endl;
		//cout << "x2 = "<<_coordinates[precedente][1]<< " y2 = " <<_coordinates[precedente][2]<< endl;
		distance2 = distance2 + pow(_coordinates[attuale][1]-_coordinates[precedente][1],2) + pow(_coordinates[attuale][2]-_coordinates[precedente][2],2);
	}
	distance2 = distance2 + pow(_coordinates[sequence.back()-1][1]-_coordinates[sequence[0]-1][1],2) + pow(_coordinates[sequence.back()-1][2]-_coordinates[sequence[0]-1][2],2);
	
	return distance2;
}

void cities::manualInitializer(vector<vector<double>> positions){

	for(int i=0; i<positions.size(); i++){

		_coordinates[i][0] = i+1;
		_coordinates[i][1] = positions[i][0];
		_coordinates[i][2] = positions[i][1];
	}
}

double cities::generateRandom(){

	return _gen.Rannyu();

}
