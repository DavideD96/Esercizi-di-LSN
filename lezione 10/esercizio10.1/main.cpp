#include <iostream>
#include <fstream>
#include "random.h"
#include "SimulatedAnnealing.h"
#include <cmath>
#include <vector>

using namespace std;

int main(){

int ncit = 32;
double Tin = 15;
double Tfin = 0.01;
double Tstep = 0.01;
double nstep = 5000;
vector<double> allT;
int Tn = static_cast<int>((Tin-Tfin)/Tstep);

for(int i=0; i<Tn; i++){
	allT.push_back(Tin-i*Tstep);
}
	
/************************** primo test di SA *************************
vector<vector<double>> cities_ = {{1,1},{4,1},{1,5},{4,5},{3.5,3},{-1,3},{1.5,6},{3.5,6}};
ncit = 8;

SimulatedAnnealing mySimulatedAnnealing(ncit,Tin);
//mySimulatedAnnealing.randomOnACircumference(1.);
mySimulatedAnnealing.manualInitializer(cities_);
double TipicalEnergy = mySimulatedAnnealing.searchTypicalEnergy(1000);
cout << TipicalEnergy << endl;
//mySimulatedAnnealing.setTemp(TipicalEnergy*10);

for(int i=0; i<Tn; i++){
	mySimulatedAnnealing.setTemp(allT[i]);
	for(int h=0; h<nstep; h++){
		mySimulatedAnnealing.move();
	}
}

/*vector<int> optPath = mySimulatedAnnealing.getPath();

for(int i=0; i<ncit; i++){
	cout << optPath[i] << " ";
}*/
//********************************************************************/

//******************** random su una circonferenza ********************

SimulatedAnnealing mySimulatedAnnealing(ncit,Tin);
mySimulatedAnnealing.randomOnACircumference(1.);

ofstream citiescirc;
citiescirc.open("circ_lenght.txt");
vector<vector<double>> circ;
vector<int> circ_;

for(int i=0; i<Tn; i++){
	mySimulatedAnnealing.setTemp(allT[i]);
	for(int h=0; h<nstep; h++){
		mySimulatedAnnealing.move();
	}
	circ_ = mySimulatedAnnealing.getPath();
	citiescirc << i << "  " << mySimulatedAnnealing.distanceTraveled(circ_) << endl;
}


citiescirc.open("circ.txt");

circ = mySimulatedAnnealing.getCoordinates();

for(int i=0; i< ncit; i++){
	citiescirc << circ[i][1] << "  " << circ[i][2] << endl; 
}


citiescirc.close();
citiescirc.open("circ_path.txt");
vector<int> optPath = mySimulatedAnnealing.getPath();

for(int i=0; i<ncit; i++){
	citiescirc << optPath[i] << endl;
}

citiescirc.close();

//********************************************************************/

//************************ random in un quadrato **********************
allT.clear();

ncit = 32;
Tin = 20;
Tfin = 0.01;
Tstep = 0.01;
nstep = 5000;

Tn = static_cast<int>((Tin-Tfin)/Tstep);

ofstream citiessquare;

SimulatedAnnealing mySimulatedAnnealing1(ncit,Tin);
mySimulatedAnnealing1.randomInASquare(1.);

citiessquare.open("square_lenght.txt");
vector<vector<double>> square;
vector<int> square_;

for(int i=0; i<Tn; i++){
	mySimulatedAnnealing1.setTemp(allT[i]);
	for(int h=0; h<nstep; h++){
		mySimulatedAnnealing1.move();
	}
	square_ = mySimulatedAnnealing1.getPath();
	citiessquare << i << "  " << mySimulatedAnnealing1.distanceTraveled(square_) << endl;
}

citiessquare.close();

citiessquare.open("square.txt");

square = mySimulatedAnnealing1.getCoordinates();

for(int i=0; i< ncit; i++){
	citiessquare << square[i][1] << "  " << square[i][2] << endl; 
}


citiessquare.close();
citiessquare.open("square_path.txt");
optPath.clear();
optPath = mySimulatedAnnealing1.getPath();

for(int i=0; i<ncit; i++){
	citiessquare << optPath[i] << endl;
}

citiessquare.close();

//*****************************************************************

return 0;
}
