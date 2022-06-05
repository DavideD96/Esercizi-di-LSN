#include "SimulatedAnnealing.h"

SimulatedAnnealing::SimulatedAnnealing(unsigned int ncit, double Tin): cities(ncit){

	_ncities = ncit;
	_temp = Tin;
	_path.resize(ncit,0);
	 vector<vector<int>> appoggio = InitializePopulation(1,ncit);
	_path = appoggio[0];
	_energy = distanceTraveled(_path);
}

double SimulatedAnnealing::searchTypicalEnergy(int nmosse){

	double sumEnergy = 0.;
	vector<vector<int>> pop = InitializePopulation(nmosse,_ncities);
	
	for(int i=0; i<nmosse; i++){
		sumEnergy += distanceTraveled(pop[i]); 
	}

	return sumEnergy/(double)nmosse;
}

void SimulatedAnnealing::move(){

	double rnd = generateRandom();

	vector<int> new_path = _path;
	//******************* permutazione ********************

	if(rnd < 0.8){

		double probability;
		double sorted;
		int ind1, ind2, appoggio;	
	
		sorted = generateRandom();
			
		if(sorted < probability){
			
			//cout << "indice permutazione = " << i << endl;
			ind1 = floor(1+(_ncities-1)*generateRandom());
			ind2 = floor(1+(_ncities-1)*generateRandom());
			appoggio = new_path[ind2];
			new_path[ind2] = new_path[ind1];
			new_path[ind1] = appoggio;
		}
	
	}else if(0.8 < rnd && rnd < 0.9){
	
	//******************* inversione ***********************

		double probability;
		double sorted;
		bool avvenuta_inversione = false;
		int m;

		int ind;
		vector<int> appoggio;

		do{
			m = floor(2+(_ncities-2)*generateRandom());
			probability = (1./static_cast<double>(m));
			sorted = generateRandom();
			
			if(sorted < probability){
	
				avvenuta_inversione = true;
	
				//cout << "indice inversione = " << i << endl;
				ind = floor(1+(_ncities-m)*generateRandom());
	
				for(int j=0; j<m; j++){
					appoggio.push_back(new_path[ind+j]);
				}
	
				for(int j=0; j<m; j++){
					new_path[ind+m-1-j] = appoggio[j];
				}
	
				appoggio.clear();	
			}
		}while(avvenuta_inversione == false);
        
	}else{

	//********************* scambio ************************
	
		double probability;
		double sorted;
		int ind1, ind2, ind_min, ind_max;
		int appoggio;
		bool avvenuto_scambio = false;
		int m;

		do{

			do{
				m = floor(2+(_ncities-2)*generateRandom());		
			}while(m > (_ncities-1)/2);
	
			probability = 1./static_cast<double>(m);
			sorted = generateRandom();
			
			if(sorted < probability){
			
					avvenuto_scambio = true;
				do{
					ind1 = floor(1+(_ncities-m)*generateRandom());
					ind2 = floor(1+(_ncities-m)*generateRandom());
	
					if(ind1 < ind2){
						ind_max = ind2;
						ind_min = ind1;
					}else{
						ind_max = ind1;
						ind_min = ind2;	
					}
			
				}while(ind_max-ind_min < m || _ncities-ind_max < m);
	
				for(int j=0; j<m; j++){
					appoggio = new_path[ind1+j];
					new_path[ind1+j] = new_path[ind2+j];
					new_path[ind2+j] = appoggio;
				}
			}
		}while(avvenuto_scambio == false);
	}
	//_path = {1,2,5,4,8,7,3,6};
	//new_path = {1,2,4,5,8,7,3,6};

	double new_energy = distanceTraveled(new_path);

	//cout << endl;
	//cout << "ne " << new_energy << endl;
	//cout << "oe " << _energy << endl;
	
	double old_p = exp(-_energy/_temp);
	double new_p = exp(-new_energy/_temp);

	double ar_probability;

	if (new_energy > _energy){
		//cout << "energia aumenta" << endl;
		ar_probability = new_p/old_p;
	}
	else{
		//cout << "energia diminuisce "<< endl;
		ar_probability = 1.;
	}

	double ar = generateRandom();
	
	if (ar < ar_probability){
		_path = new_path;
		_energy = new_energy;
	}
	/*for(int i=0; i<8; i++){
		cout << _path[i] << endl;
	}
	cout << endl;*/
}

void SimulatedAnnealing::setTemp(double T){

	_temp = T;
}

vector<int> SimulatedAnnealing::getPath(){

	return _path;
}
