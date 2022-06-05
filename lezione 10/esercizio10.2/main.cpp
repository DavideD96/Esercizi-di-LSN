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
#include "mpi.h"

using namespace std;

int main(int argc, char* argv[]){

int nprocessi = 4;
int size, rank;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Status stat1, stat2, stat3, stat4;

//************************ algoritmo completo: città su un cerchio ************************
int ncitta = 32;
int npop = 100;
int run_dim = 100000;
int run_migr = 1000;
int n_migr = run_dim/run_migr;
vector<vector<double>> circ(ncitta,vector<double>(3));
vector<double> circx(ncitta);
vector<double> circy(ncitta);

ofstream citiescirc;
//cities myCit(ncitta);

if(rank == 0){ //il nodo 0 sceglie le posizioni delle città

	cities myCit(ncitta);
	myCit.randomInASquare(1.);
	circ = myCit.getCoordinates();

	for(int i=0;i<ncitta;i++){
		circx[i] = circ[i][1];
		circy[i] = circ[i][2];
	}

	citiescirc.open("coordinates.txt");

	for(int i=0; i< ncitta; i++){
		citiescirc << circ[i][1] << "  " << circ[i][2] << endl; 
	}


	citiescirc.close();

	/*citiescirc.open("prova_dist"+to_string(rank)+".txt");
	vector<int> pathprova = {1,2,3,4,5};
	citiescirc << myCit.distanceTraveled(pathprova) << endl; 
	citiescirc.close();*/
}

int itag[nprocessi]; 

for(int i=0; i<nprocessi; i++){
	itag[i]=i+1;
}

//comunico le posizione delle città a tutti
MPI_Bcast(circx.data(), circx.size(), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
MPI_Bcast(circy.data(), circy.size(), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);

for(int i=0;i<ncitta;i++){
	circ[i][0] = i+1;
	circ[i][1] = circx[i];
	circ[i][2] = circy[i];
}

geneticAlgorithm myGA(ncitta, npop);

vector<vector<double>> circ_(ncitta, vector<double>(3));
for(int i=0; i<ncitta; i++){

	circ_[i][0] = circ[i][1];
	circ_[i][1] = circ[i][2];
}

myGA.manualInitializer(circ_);
myGA.initializePop(rank);

vector<int> bestTr(ncitta);
vector<int> bestTrApp(ncitta);
vector<int> nodi_per_scambio_(nprocessi);
vector<vector<int>> nodi_per_scambio(nprocessi);
//vector<vector<int>> pop;

for(int i=0; i<run_migr; i++){

	//scambio i percorsi migliori

	//estraggo la combinazione per gli scambi...
	if(rank==0){
		nodi_per_scambio = InitializePopulation(1,nprocessi,0);
		nodi_per_scambio_ = nodi_per_scambio[0];
	}

	//...la comunico a tutti i nodi
	MPI_Bcast(nodi_per_scambio_.data(), nodi_per_scambio_.size(), MPI_INTEGER, 0, MPI_COMM_WORLD);

	//prendo il percorso migliore per ogni continente
	bestTr = myGA.bestTravel();
	bestTrApp = bestTr;
	
	//scambio!
	if(rank == nodi_per_scambio_[0]-1){ //************************
		MPI_Send(bestTr.data(), bestTr.size(), MPI_INTEGER, nodi_per_scambio_[1]-1, itag[0], MPI_COMM_WORLD);
		MPI_Recv(bestTr.data(), bestTr.size(), MPI_INTEGER, nodi_per_scambio_[1]-1, itag[1], MPI_COMM_WORLD, &stat1);
	}else if(rank == nodi_per_scambio_[1]-1){
		MPI_Recv(bestTr.data(), bestTr.size(), MPI_INTEGER, nodi_per_scambio_[0]-1, itag[0], MPI_COMM_WORLD, &stat2);
		MPI_Send(bestTrApp.data(), bestTrApp.size(), MPI_INTEGER, nodi_per_scambio_[0]-1, itag[1], MPI_COMM_WORLD);	
	}else if(rank == nodi_per_scambio_[2]-1){
		MPI_Send(bestTr.data(), bestTr.size(), MPI_INTEGER, nodi_per_scambio_[3]-1, itag[2], MPI_COMM_WORLD);
		MPI_Recv(bestTr.data(), bestTr.size(), MPI_INTEGER, nodi_per_scambio_[3]-1, itag[3], MPI_COMM_WORLD, &stat3);
	}else if(rank == nodi_per_scambio_[3]-1){
		MPI_Recv(bestTr.data(), bestTr.size(), MPI_INTEGER, nodi_per_scambio_[2]-1, itag[2], MPI_COMM_WORLD, &stat4);
		MPI_Send(bestTrApp.data(), bestTrApp.size(), MPI_INTEGER, nodi_per_scambio_[2]-1, itag[3], MPI_COMM_WORLD);
	}
	//sostituisco nella popolazione il nuovo percorso
	myGA.replacePath(bestTr);//****************************

	/*if(rank == 0){
		for(int i=0; i<npop; i++){
			for(int j=0; j<ncitta; j++){
				myGA.pointReplace(j+1,i,j);
			}
		}
		citiescirc.open("distanze_min"+to_string(rank)+".txt",ios::app);
		pop = myGA.getPopulation();
		/*for(int i=0; i<npop; i++){
			for(int j=0; j<ncitta; j++){
				citiescirc << pop[i][j] << "  ";
			}
			citiescirc << endl;
		}*/
		/*for(int j=0; j<ncitta; j++){
				pop[0][j] = j+1;
				//citiescirc << pop[0][j] << "  ";
		}*/
		/*citiescirc << myGA.distanceTraveled(pop[0]) << endl;
		pop.clear();
		citiescirc.close();
	}*/
	//evoluzione
	for(int k=0; k<n_migr; k++){

		citiescirc.open("distanze_min"+to_string(rank)+".txt",ios::app);
		citiescirc << myGA.minLenght() << endl;
		citiescirc.close();

		myGA.selection();
		myGA.crossover();
		myGA.mutation_permutation();
		myGA.mutation_inversion_mcities();
		myGA.mutation_exchange();
	}
}

bestTr = myGA.bestTravel();
citiescirc.open("percorso_finale"+to_string(rank)+".txt");

for(int i=0; i<ncitta; i++){
	citiescirc << bestTr[i] << endl; 
}

citiescirc.close();

MPI_Finalize();

return 0;
}

