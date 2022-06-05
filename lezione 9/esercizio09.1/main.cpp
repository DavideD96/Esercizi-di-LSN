#include <iostream>
#include <fstream>
#include "random.h"
#include "checkPath.h"
#include "InitializePopulation.h"
#include "cities.h"
#include <cmath>
#include <vector>
#include <iomanip>
#include <string>
#include "geneticAlgorithm.h"

using namespace std;

int main(){

/********************************** test del checkPath ***********************************

vector<vector<int>> cammini {{1,2,3,4,5},{1,3,2,1,4},{1,2,3,4,3},{1,1,3,4,5},{1,2,3,1,1}};
bool check = false;
for(int i=0; i<5; i++){
	check = checkPath(cammini[i]);
	cout << check << " " << endl;
}

//*****************************************************************************************/

/********************************* test del generatore di cammini **************************
int ncitta = 12;
int npop = 10;

vector<vector<int>> pop = InitializePopulation(npop,ncitta);

for(int i=0; i<npop; i++){
	for(int j=0; j<ncitta; j++){
		cout << pop[i][j] << " ";	
	}
	cout << endl;
}
//*****************************************************************************************/

/************************************* test cities ****************************************
vector<vector<double>> cities_ = {{1,1},{4,1},{1,5}};
vector<int> seq = {1,3,2};

cities mycit(3);
mycit.manualInitializer(cities_);
double dist = mycit.distanceTraveled(seq);

cout << dist << endl;
//*****************************************************************************************/

/************************************* test selezione *************************************
int ncitta = 5;
int npop = 7;

vector<vector<double>> cities_ = {{1,1},{4,1},{1,5},{4,5},{3,3}};

geneticAlgorithm myGA(ncitta, npop);
myGA.manualInitializer(cities_);
myGA.initializePop();

vector<vector<int>> pop = myGA.getPopulation();

for(int i=0; i<npop; i++){
	for(int k=0; k<ncitta; k++){
		cout << pop[i][k] << " ";
	}
	cout << endl;
}

myGA.selection();
pop = myGA.getPopulation();

cout << " selezionata: " << endl;
for(int i=0; i<npop; i++){
	for(int k=0; k<ncitta; k++){
		cout << pop[i][k] << " ";
	}
	cout << endl;
}
//****************************************************************************************/

/************************************* test mutazioni ************************************
int ncitta = 8;

int npop = 100;

vector<vector<double>> cities_ = {{1,1},{4,1},{1,5},{4,5},{3,3},{-1,3},{1.5,6},{3.5,6}};


geneticAlgorithm myGA(ncitta, npop);
myGA.manualInitializer(cities_);
myGA.initializePop();


vector<vector<int>> pop = myGA.getPopulation();

for(int i=0; i<npop; i++){

	cout << i << "   ";

	for(int k=0; k<ncitta; k++){
		cout << pop[i][k] << " ";
	}
	cout << endl;
}


myGA.mutation_permutation();
myGA.mutation_inversion_mcities();
pop = myGA.getPopulation();

cout << " selezionata: " << endl;

for(int i=0; i<npop; i++){

	cout << i << "   ";

	for(int k=0; k<ncitta; k++){
		cout << pop[i][k] << " ";
	}

	cout << endl;
}
//********************************************************************************************/

/*************************** la mia prima ricerca! (solo mutazioni) **************************
int ncitta = 8;
int npop = 100;
int run_dim = 1000;

vector<vector<double>> cities_ = {{1,1},{4,1},{1,5},{4,5},{3,3},{-1,3},{1.5,6},{3.5,6}};

geneticAlgorithm myGA(ncitta, npop);
myGA.manualInitializer(cities_);
myGA.initializePop();

vector<vector<int>> pop = myGA.getPopulation();

for(int k=0; k<run_dim; k++){

	myGA.selection();
	myGA.mutation_permutation();
	myGA.mutation_inversion_mcities();
}

pop = myGA.getPopulation();

cout << " ultima generazione: " << endl;
for(int i=0; i<npop; i++){

	cout << i << "   ";

	for(int k=0; k<ncitta; k++){
		cout << pop[i][k] << " ";
	}
	cout << endl;
}
//****************************************************************************************/

/*********************************** test del crossover **********************************

int ncitta = 8;
int npop = 20;

vector<vector<double>> cities_ = {{1,1},{4,1},{1,5},{4,5},{3,3},{-1,3},{1.5,6},{3.5,6}};


geneticAlgorithm myGA(ncitta, npop);
myGA.manualInitializer(cities_);
myGA.initializePop();


vector<vector<int>> pop = myGA.getPopulation();

for(int i=0; i<npop; i++){


	cout << i << "   ";


	for(int k=0; k<ncitta; k++){
		cout << pop[i][k] << " ";
	}

	cout << endl;
}

//myGA.selection();
myGA.crossover();

pop = myGA.getPopulation();


cout << " dopo il crossover: " << endl;
for(int i=0; i<npop; i++){

	cout << i << "   ";


	for(int k=0; k<ncitta; k++){
		cout << pop[i][k] << " ";
	}

	cout << endl;
}
//****************************************************************************************/

/************************************* test completo *************************************
int ncitta = 8;
int npop = 100;
int run_dim = 1000;

vector<vector<double>> cities_ = {{1,1},{4,1},{1,5},{4,5},{3,3},{-1,3},{1.5,6},{3.5,6}};


geneticAlgorithm myGA(ncitta, npop);
myGA.manualInitializer(cities_);
myGA.initializePop();


vector<vector<int>> pop = myGA.getPopulation();

for(int k=0; k<run_dim; k++){

	myGA.selection();
	myGA.crossover();
	myGA.mutation_permutation();
	myGA.mutation_inversion_mcities();
	myGA.mutation_exchange();

}

pop = myGA.getPopulation();


cout << " ultima generazione: " << endl;
for(int i=0; i<npop; i++){

	cout << i << "   ";


	for(int k=0; k<ncitta; k++){
		cout << pop[i][k] << " ";
	}

	cout << endl;
}
//*****************************************************************************************/

//************************ algoritmo completo: città su un cerchio ************************
int ncitta = 32;
int npop = 100;
int run_dim = 2000;

ofstream citiescirc;
citiescirc.open("circ.txt");

geneticAlgorithm myGA(ncitta, npop);
myGA.randomOnACircumference(1.);
myGA.initializePop();

vector<vector<double>> circ;

circ = myGA.getCoordinates();

for(int i=0; i< ncitta; i++){
	citiescirc << circ[i][1] << "  " << circ[i][2] << endl; 
}

citiescirc.close();

//vector<vector<int>> pop = myGA.getPopulation();
/*cout << "popolazione iniziale: " << endl;

for(int i=0; i<npop; i++){

	cout << i << "   ";

	for(int k=0; k<ncitta; k++){
		cout << pop[i][k] << " ";
	}
	cout << endl;
}*/

citiescirc.open("circ_meanL.txt");
double meanLenght, minLenght;

for(int k=0; k<run_dim; k++){

	meanLenght = myGA.meanBestHalf();
	minLenght = myGA.minLenght();
	citiescirc << k << "  " << meanLenght << "  " << minLenght << endl;

	myGA.selection();
	myGA.crossover();
	myGA.mutation_permutation();
	myGA.mutation_inversion_mcities();
	myGA.mutation_exchange();

}

citiescirc.close();
//pop = myGA.getPopulation();

/*cout << endl << "popolazione finale: " << endl;
for(int i=0; i<npop; i++){

	cout << i << "   ";

	for(int k=0; k<ncitta; k++){
		cout << pop[i][k] << " ";
	}
	cout << endl;
}*/

vector<int> bestTr = myGA.bestTravel();

citiescirc.open("circ_path.txt");
cout << "Best travel" << endl;
for(int k=0; k<ncitta; k++){
	cout << bestTr[k] << " ";
	citiescirc << bestTr[k] << endl;
}
cout << endl;
citiescirc.close();
//****************************************************************************************/

//*********************** algoritmo completo: città in un quadrato ************************
ncitta = 32;
npop = 100;
run_dim = 20000;

geneticAlgorithm myGA1(ncitta, npop);
myGA1.randomInASquare(1.);
myGA1.initializePop();

ofstream citiessquare;
citiessquare.open("square.txt");

vector<vector<double>> square;

square = myGA1.getCoordinates();

for(int i=0; i< ncitta; i++){
	citiessquare << square[i][1] << "  " << square[i][2] << endl; 
}

citiessquare.close();
//vector<vector<int>> pop = myGA1.getPopulation();
//cout << "popolazione iniziale: " << endl;


/*for(int i=0; i<npop; i++){

	cout << i << "   ";


	for(int k=0; k<ncitta; k++){
		cout << pop[i][k] << " ";
	}
	cout << endl;

}*/
citiessquare.open("square_meanL.txt");

for(int k=0; k<run_dim; k++){

	meanLenght = myGA1.meanBestHalf();
	minLenght = myGA1.minLenght();
	citiessquare << k << "  " << meanLenght << "  " << minLenght << endl;

	myGA1.selection();
	myGA1.crossover();
	myGA1.mutation_permutation();
	myGA1.mutation_inversion_mcities();

	myGA1.mutation_exchange();

}

citiessquare.close();
//pop = myGA1.getPopulation();


//cout << endl << "popolazione finale: " << endl;
/*for(int i=0; i<npop; i++){


	cout << i << "   ";

	for(int k=0; k<ncitta; k++){
		cout << pop[i][k] << " ";

	}
	cout << endl;
}*/


vector<int> bestTr1 = myGA1.bestTravel();

citiessquare.open("square_path.txt");
cout << "Best travel" << endl;
for(int k=0; k<ncitta; k++){

	cout << bestTr1[k] << " ";
	citiessquare << bestTr1[k] << endl;
}
cout << endl;
citiessquare.close();
//****************************************************************************************/

return 0;
}

