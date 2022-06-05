#include "geneticAlgorithm.h"

geneticAlgorithm::geneticAlgorithm(unsigned int ncit, unsigned int npop): cities(ncit){

	_ncities = ncit;
	_npopulation = npop;

	_singlePath.resize(ncit,0);
	_distances.resize(npop,0);
	_population.resize(npop,_singlePath);
	
}

void geneticAlgorithm::initializePop(){

	_population = InitializePopulation(_npopulation, _ncities);
	
}

double geneticAlgorithm::meanBestHalf(){

	//*************** ordinamento ****************
	for(int i=0; i<_npopulation; i++){
		_distances[i] = distanceTraveled(_population[i]);	
	}
	
	vector<double> distancesAus = _distances;
	vector<double> distancesInd(_npopulation);
	vector<vector<int>> orderedPopulation;
	int max = 10000;
	int min = max;
	int ind = 0;
	for(int i=0; i<_npopulation; i++){

		min = max;
		for (int j=0; j<_npopulation; j++){

			if(min > distancesAus[j]){
				min = distancesAus[j];
				ind = j;
			}
		}
		distancesAus[ind] = max+1;
		distancesInd[i] = ind;
		orderedPopulation.push_back(_population[ind]);
	}	
	vector<vector<int>> bestHalf(_ncities/2 + 1);

	auto in = orderedPopulation.begin();
	auto fin = orderedPopulation.begin() + _ncities/2 + 1;
	
	copy(in, fin, bestHalf.begin());
	double dist = 0.;	

	for(int i=0; i<_ncities/2; i++){
		dist += distanceTraveled(bestHalf[i]);
	}

	return dist*2./_ncities;
}

double geneticAlgorithm::minLenght(){

	vector<int> bestTr = bestTravel();
	return distanceTraveled(bestTr);
}

void geneticAlgorithm::selection(){


	/*for(int i=0; i<_npopulation; i++){
		for(int k=0; k<_ncities; k++){
			cout << _population[i][k] << " ";
		}
		cout << endl;
	}*/

	//*************** ordinamento ****************
	for(int i=0; i<_npopulation; i++){
		/*for(int j=0; j<_ncities; j++){
			cout << _population[i][j] << " ";
		}
		cout << endl;*/
		_distances[i] = distanceTraveled(_population[i]);	
	}
	//cout << endl;
	vector<double> distancesAus = _distances;
	vector<double> distancesInd(_npopulation);
	vector<vector<int>> orderedPopulation;//(_npopulation, _ncities);
	int max = 10000;
	int min = max;
	int ind = 0;
	for(int i=0; i<_npopulation; i++){

		min = max;
		for (int j=0; j<_npopulation; j++){

			if(min > distancesAus[j]){
				min = distancesAus[j];
				ind = j;
			}
		}
		distancesAus[ind] = max+1;
		distancesInd[i] = ind;
		orderedPopulation.push_back(_population[ind]);
		
		/*for(int k=0; k<_ncities; k++){
			cout << orderedPopulation[i][k] << " ";
		}
		cout << endl;*/
	}
	//*****************************************************
	
	//******************** selezione ********************** 
	double r;
	double p = 1.4;
	int indice;

	for(int i=0; i<_npopulation; i++){

		r = generateRandom();
		indice = floor(_npopulation*pow(r,p));
		//cout << indice << endl;
		_population[i] = orderedPopulation[indice];
		//cout << indice << endl;
		
	}
	//*****************************************************
	
}

vector<vector<int>> geneticAlgorithm::getPopulation(){

	return _population;

}

void geneticAlgorithm::mutation_permutation(){
	
	double probability = 0.08;
	double sorted;
	int ind1, ind2, appoggio;	

	for(int i=0; i<_npopulation; i++){

		sorted = generateRandom();
		
		if(sorted < probability){
			
			//cout << "indice permutazione = " << i << endl;
			ind1 = floor(1+(_ncities-1)*generateRandom());
			ind2 = floor(1+(_ncities-1)*generateRandom());
			appoggio = _population[i][ind2];
			_population[i][ind2] = _population[i][ind1];
			_population[i][ind1] = appoggio;
		}
	}
}

void geneticAlgorithm::mutation_inversion_mcities(){

	double probability;
	double sorted;
	int m;

	//cout << "inizio inversion ***************************** " << endl;
	int ind;
	vector<int> appoggio;

	for(int i=0; i<_npopulation; i++){

		m = floor(2+(_ncities-2)*generateRandom());
		probability = 0.14*(1./static_cast<double>(m));
		sorted = generateRandom();
		
		if(sorted < probability){

			//cout << "indice inversione = " << i << endl;
			ind = floor(1+(_ncities-m)*generateRandom());

			for(int j=0; j<m; j++){
				appoggio.push_back(_population[i][ind+j]);
			}

			for(int j=0; j<m; j++){
				_population[i][ind+m-1-j] = appoggio[j];
			}

			appoggio.clear();	
		}
	}
	//cout << "fine inversion ***************************** " << endl;
}

void geneticAlgorithm::mutation_exchange(){

	//cout << "inizio exchange ***************************** " << endl;
	double probability;
	double sorted;
	int ind1, ind2, ind_min, ind_max;
	int appoggio;
	int m;

	for(int i=0; i<_npopulation; i++){

		do{
			m = floor(2+(_ncities-2)*generateRandom());		
		}while(m > (_ncities-1)/2);

		probability = 0.14*(1./static_cast<double>(m));
		sorted = generateRandom();
		
		if(sorted < probability){

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
				appoggio = _population[i][ind1+j];
				_population[i][ind1+j] = _population[i][ind2+j];
				_population[i][ind2+j] = appoggio;
			}
		}
	}
	//cout << "fine exchange ***************************** " << endl;
}

void geneticAlgorithm::crossover(){

	double probability = 0.50;
	double sorted;
	int ind;
	int min;
	int path1 = 1;
	int path2 = 1;
	int npop_disponibile = _npopulation;
	vector<int> indices1, indices2;
	vector<int> appoggio1, appoggio2;
	vector<vector<int>> newgeneration;// = _population;

	//cout << "population in crossover" << endl;

	/*for(int l=0; l<_npopulation; l++){
		for(int k=0; k<_ncities; k++){
			cout << newgeneration[l][k] << " ";
		}
		cout << endl;
	}*/

	for(int i=0; i<_npopulation/2; i++){ //nota la divisione fra interi!

		sorted = generateRandom();
		
		if(sorted < probability){

			do{
				path1 = floor((npop_disponibile)*generateRandom());
				path2 = floor((npop_disponibile)*generateRandom());

			}while(path1 == path2);

			//cout << "crossover tra: " << path1 << " " << path2 << endl;

			ind = floor(1+(_ncities-2)*generateRandom());

			//cout << "indice = " << ind << endl;
			indices1.resize(_ncities-ind,0);
			indices2.resize(_ncities-ind,0);
			//cout << "path2 " << endl;

			for(int j=ind; j<_ncities; j++){ 
			
				//modifico path2 sulla base di path1
				indices1[j-ind] = find(_population[path1],_population[path2][j]);

				//cout << _population[path2][j] << " ";

				//modifico path1 sulla base di path2
				indices2[j-ind] = find(_population[path2],_population[path1][j]);
			}
			//cout << endl;
			//cout << endl;
			/*cout << "indices" << endl;
			for(int j=0; j<_ncities-ind; j++){ 
			
				cout << indices1[j] << " ";
			}*/			

			for(int k=0; k<ind; k++){
				
				appoggio1.push_back(_population[path1][k]);
				appoggio2.push_back(_population[path2][k]);
			}

			for(int j=0; j<_ncities-ind; j++){ 
			
				min = findMin(indices1);
				//cout << "new gen "<< path2<< " indice " <<ind+j<< " "<< newgeneration[path2][ind+j] << endl;
				//cout << "sostituto "<<_population[path1][indices1[min]]<<endl;
				appoggio2.push_back(_population[path1][indices1[min]]);
				indices1.erase(indices1.begin()+min);

				min = findMin(indices2);
				appoggio1.push_back(_population[path2][indices2[min]]);
				indices2.erase(indices2.begin()+min);
			}			

			newgeneration.push_back(appoggio1);
			newgeneration.push_back(appoggio2);

			appoggio1.clear();
			appoggio2.clear();
			/*for(int l=0; l<npop_disponibile; l++){
				cout << l << "   ";
				for(int k=0; k<_ncities; k++){
					cout << _population[l][k] << " ";
				}
				cout << endl;
			}*/

			//cout << "vecchio " << endl;
			/*for(int l=0; l<npop_disponibile; l++){
				cout << l << "   ";
				for(int k=0; k<_ncities; k++){
					cout << _population[l][k] << " ";
				}
				cout << endl;
			}*/
			//cout << "nuovo " << endl;
			/*for(int l=0; l<newgeneration.size(); l++){
				cout << l << "   ";
				for(int k=0; k<_ncities; k++){
					cout << newgeneration[l][k] << " ";
				}
				cout << endl;
			}*/

			//cout << endl;	

			_population.erase(_population.begin()+path1);

			if(path1 < path2){
				_population.erase(_population.begin()+path2-1);
			}else{
				_population.erase(_population.begin()+path2);
			}

			npop_disponibile = npop_disponibile - 2;
		}	
	}
	
	for(int i=0; i<_population.size(); i++){
		newgeneration.push_back(_population[i]); //aggiungo i percorsi che non hanno fatto crossover
	}

	_population = newgeneration;
}

vector<int> geneticAlgorithm::bestTravel(){
	
	double minDistance = distanceTraveled(_population[0]);
	int index = 0;

	for(int i=1; i<_npopulation; i++){
		
		if(distanceTraveled(_population[i])<minDistance){

			minDistance = distanceTraveled(_population[i]);
			index = i;
		}
	}
	return _population[index];
}

int geneticAlgorithm::find(vector<int> vett,int target){

	int ind;
	//cout << "******* find ********" << endl;

	for(int i=0; i<vett.size(); i++){

		//cout << vett[i] << endl;
		if(vett[i] == target){

			ind = i;
		}
	}
	//cout << "ind = " << ind << endl;
	//cout << "*********************" << endl;
	return ind;
}

int geneticAlgorithm::findMin(vector<int> vett){

	int min = 0;
	//cout << "******* find min ********" << endl;


	for(int i=0; i<vett.size(); i++){
		//cout << vett[i] << endl;
		if(vett[min]>vett[i]){
			min = i;
		}
	}
	//cout << "min = " << min << endl;
	//cout << "*********************" << endl;
	return min;
}
