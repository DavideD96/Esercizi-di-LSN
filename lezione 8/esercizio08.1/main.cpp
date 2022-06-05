#include <iostream>
#include <fstream>
#include "random.h"
#include "Metropolis.h"
#include "probDistr.h"
#include "probDensityTrial.h"
#include <cmath>
#include <vector>

using namespace std;

double Hpsi_over_psi(double, double, double);
double error(unsigned int, double, double);

int main(){

int distribution_type = 0; //0 means T is an uniform distribution, 1 means T is a gaussian distribution.
double mu = 1.;
double sigma = 1.;

probDensityTrial* myPDT = new probDensityTrial();

//Initial position: centre of the atom (is surely the best choice).
Metropolis myMetropolis(myPDT, 0, distribution_type);

//set parameters
myPDT->setParam(mu,sigma);

//Equilibration.
for(int i=0; i<1000; i++){
	myMetropolis.move();
}

double position;

int nBlocks = 100;
int nInBlocks = 10000;
double mean_Energy;
double mean_Energy2;
double sum_Energy = 0;

double AccAtt = 0.;

double sum = 0.;
double sum2 = 0.;

vector<double> estimations_Energy(nBlocks,0);
vector<double> estimations_uncertanties(nBlocks,0);

for(int j=0; j<nBlocks; j++){

	sum_Energy = 0;

	for(int i=0; i<nInBlocks; i++){
		myMetropolis.move();
		position = myMetropolis.get_position();
		sum_Energy += Hpsi_over_psi(position, mu, sigma);
		//cout << " x = " <<position[0]<<"  "<<" y = " <<position[1]<<"  "<<" z = " <<position[2]<<endl;
	}
	
	mean_Energy = sum_Energy/nInBlocks;
	mean_Energy2 = pow(mean_Energy,2);

	sum += mean_Energy;
	sum2 += mean_Energy2;

	estimations_Energy[j] = sum/(j+1);
	estimations_uncertanties[j] = error(j+1, sum/(j+1), sum2/(j+1));

	AccAtt = myMetropolis.get_Acc_rate();
	myMetropolis.reset_Acc_rate();
	cout << "Acc/Att ratio = " << AccAtt << "   " << "Mean = " << estimations_Energy[j] << "  Error = " << estimations_uncertanties[j] << endl;



}

return 0;
}

double Hpsi_over_psi(double x, double mu, double sigma){

	double psi = exp(-pow(x-mu,2)/(2*sigma*sigma))+exp(-pow(x+mu,2)/(2*sigma*sigma));
	double Vpsi = (pow(x,4)-(5./2.)*pow(x,2))*psi;
	double Tpsi = -0.5*(1/(sigma*sigma))*((pow(x-mu,2)/(sigma*sigma)-1)*exp(-pow(x-mu,2)/(2*sigma*sigma))+(pow(x+mu,2)/(sigma*sigma)-1)*exp(-pow(x+mu,2)/(2*sigma*sigma))); //m, htagliato = 1
	return (Vpsi+Tpsi)/psi;
}

//funzione per calcolare l'errore statistico
double error(unsigned int N, double mean, double mean2){
	if (N==1){
		return 0;
	}
	else {
		return sqrt((1./(N-1))*(mean2-pow(mean,2)));
	}
}
