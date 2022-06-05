#include "InitializePopulation.h"

vector<vector<int>> InitializePopulation(unsigned int npop, unsigned int ncit, int rank){
	
	vector<vector<int>> population(npop);

	//********************************* inizializzo generatore *******************************
	Random rand;

	ifstream inData;

	int seed[4];
	int p[8];

	inData.open("Primes.txt");

	if (inData.is_open()){
		inData >> p[0] >> p[1] >> p[2] >> p[3] >> p[4] >> p[5] >> p[6] >> p[7];
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
			rand.SetRandom(seed, p[2*rank], p[2*rank+1]);
			}
		}
		inData.close();
	}
	

	//************************************************************************
	double nrandom;
	bool checkRepetition;
	int interval;

	for(int i=0; i<npop; i++){
		
		population[i] = vector<int>(ncit);
		population[i][0] = 1; //la città di partenza è fissata

		for(int j=1; j<ncit; j++){
			
			do{
				nrandom = rand.Rannyu();
				interval = intervalFinder(nrandom,ncit);
				checkRepetition = true;

				for(int k=0; k<j; k++){

					if(population[i][k] == interval){
						checkRepetition = false;
					}
				}

			}while(checkRepetition == false);
			
			population[i][j] = interval;
		}
	}
	return population;
}

int intervalFinder(double number, unsigned int nIntervals){

	double domain_selected[2] = {0,1};
	int _K = nIntervals;
	double interval_lenght = (domain_selected[1]-domain_selected[0])/_K;
	int range_int[2] = {1, _K}; //intervalli numerati da 1 a _K; man mano che vado avanti con la ricerca riduco questo range.
	int INT_NUMBER = 0; //numero dell'intervallo
	unsigned int k_ = _K;
	
	INT_NUMBER = 0;

	while (INT_NUMBER == 0){
		if (k_ % 2 == 0){									//PARI
			if (number < (domain_selected[1]+domain_selected[0])/2.){		//se vale questa...
				domain_selected[1] = (domain_selected[1]+domain_selected[0])/2.;	//allora devo cercare nella metà di sinistra degli intervalli!
				range_int[1] = range_int[1] - k_/2;
			} else if (number == (domain_selected[1]+domain_selected[0])/2.){
				INT_NUMBER = range_int[1] - k_/2 + 1; 					//se il numero è sul bordo appartiene all'intervallo a destra.
			} else {
				domain_selected[0] = (domain_selected[1]+domain_selected[0])/2.;
				range_int[0] = range_int[1] - k_/2 + 1;
			}
			k_ = k_/2;
		} else {										//DISPARI
			if (number < (domain_selected[1]+domain_selected[0])/2.){
				domain_selected[1] = (domain_selected[1]+domain_selected[0])/2.-interval_lenght/2.;
				range_int[1] = range_int[1] - k_/2 - 1;				//ricorda che k_/2 se k_ è dispari...	
				if (domain_selected[1] <= number){
					INT_NUMBER = range_int[1] + 1;
				}
			} else {
				domain_selected[0] = (domain_selected[1]+domain_selected[0])/2.+interval_lenght/2.;
				range_int[0] = range_int[1] - k_/2 + 1;	
				if (number < domain_selected[0]){
					INT_NUMBER = range_int[0] - 1;
				}
			}
			k_ = k_/2;	//ricorda che k_/2 se k_ è dispari...
		}
		if (range_int[0]==range_int[1] && INT_NUMBER==0){
			INT_NUMBER = range_int[0];
		}
	}
	return INT_NUMBER;
}
