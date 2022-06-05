/*
Autore: Davide Decastri
Programma per testare un generatore di numeri pseudo casuali.
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "random.h"
#include "distributionAnalyzer.h"

using namespace std;

double error(unsigned int, double, double);

int main(){

ifstream inData;
ofstream outData;

int seed[4];
int p1, p2;

inData.open("Primes.txt");

if (inData.is_open()){
	inData >> p1 >> p2;
} 
else{
	cerr << "Unable to open Primes" <<endl;
}
inData.close();

Random myRandom; //definisco il mio generatore di numeri casuali 

inData.open("seed.in");
string property;

if (inData.is_open()){
	while(!inData.eof()){
		inData >> property;
		if(property == "RANDOMSEED"){
			inData >> seed[0] >> seed[1] >> seed[2] >>seed[3];
			myRandom.SetRandom(seed, p1, p2);
		}
	}
	inData.close();
}

unsigned int M = 1000000; //lanci totali
unsigned int N = 100; //numero blocchi

double* randNumbers = new double[M];

//carico vettore con numeri casuali in [0;1)
for (int i=0; i<M; i++){
	randNumbers[i] = myRandom.Rannyu();
	//cout << randNumbers[i] <<endl;
}

//divido in blocchi
unsigned int L = M/N; //!!!! divisione fra interi = intero

//******************************************** media ********************************************

//calcolo delle incertezze sulla media
double sum_in_block;
double* averages = new double[N];
double* averages2 = new double[N];

for (int i=0; i<N; i++){
	sum_in_block = 0;
	for (int j=0; j<L; j++){
		sum_in_block+=randNumbers[i*L+j];
	}
	averages[i]=sum_in_block/L;
	//cout<<averages[i]<<endl;
	averages2[i]=pow(averages[i],2);
}

myRandom.SaveSeed(); //salvo il seed in modo da usare sempre lo stesso facendo esecuzioni successive.

//ora faccio la media sui tanti esperimenti Monte Carlo

double* partial_mean_arr = new double[N];
double* partial_mean2_arr = new double[N];
double* errors = new double[N];

double partial_sum = 0.;
double partial_sum2 = 0.;

ofstream WriteMean;
WriteMean.open("mean.txt");

for (int i=0; i<N; i++){
	partial_sum += averages[i];
	partial_sum2 += pow(averages[i],2);

	//riempio gli array
	partial_mean_arr[i] = partial_sum/(i+1); //media fatta sui risultati di i+1 blocchi
	partial_mean2_arr[i] = partial_sum2/(i+1); //media^2 fatta sui risultati di i+1 blocchi
	errors[i] = error(i+1, partial_mean_arr[i], partial_mean2_arr[i]);
	//cout<<left<<"Media ed errore sulla media a "<<i+1<<" blocchi: "<<endl;
	//cout<<right<<partial_mean_arr[i]<<"  "<<errors[i]<<endl;
	WriteMean<<i+1<<"  "<<partial_mean_arr[i]<<"  "<<errors[i]<<endl;

}
WriteMean.close();

//********************************** varianza **************************************

ofstream WriteVar;
WriteVar.open("var.txt");

//calcolo delle incertezze sulla varianza
double sum_in_block2;
double* averages_2 = new double[N];
double* averages_22 = new double[N];

for (int i=0; i<N; i++){
	sum_in_block2 = 0;
	for (int j=0; j<L; j++){
		sum_in_block2+=pow(static_cast<double>(randNumbers[i*L+j])-0.5,2);
	}
	averages_2[i]=sum_in_block2/L;
	//cout<<i<<endl;
	//cout<<averages_2[i]<<endl;
	averages_22[i]=pow(averages_2[i],2);
	//cout<<averages_22[i]<<endl<<endl;

}

double* partial_var_arr = new double[N];
double* partial_var2_arr = new double[N];
double* errors2 = new double[N];

partial_sum = 0.;
partial_sum2 = 0.;

for (int i=0; i<N; i++){
	partial_sum += averages_2[i];
	partial_sum2 += pow(averages_2[i],2);

	//riempio gli array
	partial_var_arr[i] = partial_sum/(i+1); //varianza fatta sui risultati di i+1 blocchi
	partial_var2_arr[i] = partial_sum2/(i+1); //varianza^2 fatta sui risultati di i+1 blocchi
	errors2[i] = error(i+1, partial_var_arr[i], partial_var2_arr[i]);
	//cout<<left<<"varianza ed errore sulla varianza a "<<i+1<<" blocchi: "<<endl;
	//cout<<right<<partial_var_arr[i]<<"  "<<errors2[i]<<endl;
	WriteVar<<i+1<<"  "<<partial_var_arr[i]<<"  "<<errors2[i]<<endl;
}

WriteVar.close();

//******************************* CHI QUADRO ***********************************

ofstream WriteChi;
WriteChi.open("chi.txt");

distributionAnalyzer myDistributionAnalyzer(randNumbers, M); //costruisco analizzatore

myDistributionAnalyzer.setIntervals(100); //divido in 100 intervallini
myDistributionAnalyzer.setPeriod(10000); //numero di estrazioni con cui calcolo il chi quadro volta per volta.
myDistributionAnalyzer.setExtrema(0.,1.);

myDistributionAnalyzer.dataInIntervals();
//int** conteggi_per_bin = myDistributionAnalyzer.getCountsInInt();
double* chiSquare = myDistributionAnalyzer.chi2();

for(int j=0; j<M/10000; j++){
	//cout<<j<<endl;
	//cout<<chiSquare[j]<<endl;
	//cout<<endl;
	WriteChi<<j+1<<"  "<<chiSquare[j]<<endl;
}

WriteChi.close();

delete[] chiSquare;

delete[] averages;
delete[] averages2;

delete[] averages_2;
delete[] averages_22;

delete[] partial_mean_arr;
delete[] partial_mean2_arr;

delete[] errors;
delete[] errors2;

return 0;
}

//funzione per calcolare l'errore statistico (?)
double error(unsigned int N, double mean, double mean2){
	if (N==1){
		return 0;
	}
	else {
		return sqrt((1./(N-1))*(mean2-pow(mean,2)));
	}
}
