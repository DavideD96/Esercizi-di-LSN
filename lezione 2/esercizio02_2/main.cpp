#include <iostream>
#include <fstream>
#include "cubicRW.h"
#include "isotropicRW.h"
#include <iomanip>

using namespace std;

//funzione che aggiorna la distanza percorsa in più all'aumentare delle simulazioni, 
//per ogni passo. Per conoscere <r^2> bisogna poi dividere per il numero di 
//simulazioni. Si considera che nella prima colonna di journey ci sia la posizione 
//iniziale (fissata).

template<typename T> //ASSICURATI CHE SIA GIUSTO USARE I TEMPLATE COSì, con il tipo specificato alla chiamata in <...>
void addDistances(double* dist2, T** journ, int npassi, int dim){

	double distanzaAggiunta;

	for(int i = 0; i<npassi; i++){
		distanzaAggiunta = 0.;
		for(int j = 0; j<dim; j++){
			distanzaAggiunta += pow(journ[j][i+1],2); //non considero 0 passi
		}
		dist2[i] += distanzaAggiunta;	
	}
}

int main(){


unsigned int nWalks = 10000;
unsigned int nBlocks = 100;
unsigned int nWalksInBlock = nWalks/nBlocks;
unsigned int nStep = 100; //journey includerà anche la posizione iniziale!

//**************************** RW CUBICO *************************
ofstream writeCubic;

int startingPoint[3] = {0,0,0};
int** cubicJourney; //***= new int*[3]; //NOTA IMPORTANTE: questa variabile (puntatore ad array di puntatori) diverrà membro della classe
double** distances2 = new double*[nBlocks]; //punto di partenza non incluso, se vuoi cambialo.
double* distances = new double[nStep]();
double* errors = new double[nStep]();

for (int i=0; i<nBlocks; i++){
	distances2[i] = new double[nStep]();
	//distances = new double[nStep]();
}
//***for(int i=0; i<3; i++){
//***	cubicJourney[i] = new int[nStep+1]; // //*** = non necessario!!!
//***}

cubicRW* myCubicRW = new cubicRW(startingPoint, 3, nStep);
//cout << distances2[99][100] << endl;

//simulo i random walk
for(int i=0; i<nBlocks; i++){
 	for(int j=0; j<nWalksInBlock; j++){
		myCubicRW->walk();
		cubicJourney = myCubicRW->getJourney();
		addDistances<int>(distances2[i], cubicJourney, nStep, 3);
	}
}

ofstream journ;
journ.open("cubicJourn.dat");

for(int i=0; i<nStep; i++){
	journ << setw(20) << cubicJourney[0][i] << setw(20) << cubicJourney[1][i] << setw(20) << cubicJourney[2][i] << endl;
}

journ.close();

double quantity = 0;
double quantity2 = 0;

//calcolo valor medio ed errore della distanza, per ogni step.
for(int i=0; i<nStep; i++){

	quantity = 0;
	quantity2 = 0;

	for(int j=0; j<nBlocks; j++){

		quantity2 += distances2[j][i]/nWalksInBlock;
		quantity += sqrt(distances2[j][i]/nWalksInBlock);

	}

	distances[i] = quantity/nBlocks;
	errors[i] = sqrt((1./(nBlocks-1))*((1./nBlocks)*quantity2-pow((1./nBlocks)*quantity,2)));
}

writeCubic.open("cubic.txt");

for(int j=0; j<nStep; j++){
	writeCubic<<j+1<<" "<<distances[j]<<" "<<errors[j]<<endl;
}

writeCubic.close();

//cout<<endl<<endl;

//************************* RW ISOTROPO *************************
ofstream writeIso;

double startingPoint_[3] = {0.,0.,0.};
double** isoJourney; //*** = new double*[3]; //NOTA IMPORTANTE: questa variabile (puntatore ad array di puntatori) diverrà membro della classe
double** distances2_ = new double*[nBlocks]; //punto di partenza non incluso, se vuoi cambialo.
double* distances_ = new double[nStep]();  //uso array diversi da quelli di prima per mantenere l'info sulla prima parte
double* errors_ = new double[nStep]();

for (int i=0; i<nBlocks; i++){
	distances2_[i] = new double[nStep]();
	//distances = new double[nStep]();
}

//***for(int i=0; i<3; i++){
//***	isoJourney[i] = new double[nStep+1];
//***}
//cout << "ciao" << endl;
isotropicRW* myIsoRW = new isotropicRW(startingPoint_, nStep);

//simulo i random walk
for(int i=0; i<nBlocks; i++){
	//cout << i << endl;
	for(int j=0; j<nWalksInBlock; j++){
		myIsoRW->walk();
		isoJourney = myIsoRW->getJourney();
		addDistances<double>(distances2_[i], isoJourney, nStep, 3);
	}
}
//cout << "ciao" << endl;
journ.open("isoJourn.dat");

for(int i=0; i<nStep; i++){
	journ << setw(20) << isoJourney[0][i] << setw(20) << isoJourney[1][i] << setw(20) << isoJourney[2][i] << endl;
}

journ.close();

//calcolo valor medio ed errore della distanza, per ogni step.
for(int i=0; i<nStep; i++){

	quantity = 0;
	quantity2 = 0;

	for(int j=0; j<nBlocks; j++){

		quantity2 += distances2_[j][i]/nWalksInBlock;
		quantity += sqrt(distances2_[j][i]/nWalksInBlock);

	}

	distances_[i] = quantity/nBlocks;
	errors_[i] = sqrt((1./(nBlocks-1))*((1./nBlocks)*quantity2-pow((1./nBlocks)*quantity,2)));
}

writeIso.open("isotropic.txt");

for(int j=0; j<nStep; j++){
	writeIso<<j+1<<" "<<distances_[j]<<" "<<errors_[j]<<endl;
}

writeIso.close();

//calcolo la radice della media della distanza percorsa per ogni step.
//for(int j=0; j<nStep; j++){
//	distances_[j] = sqrt(distances2_[j]/nWalks);
//	cout<<"Step number "<<j+1<<": "<<distances_[j]<<endl;
//}

delete myCubicRW;
delete myIsoRW;
for(int i=0; i<nBlocks; i++){
	delete[] distances2[i];
	delete[] distances2_[i];
}
delete[] distances2;
delete[] distances;
delete[] errors;
delete[] errors_;
delete[] distances2_;
delete[] distances_;

//NOTA cubicJourney e isoJourney vengono distrutti attraverso il distruttore della classe.

return 0;
}


