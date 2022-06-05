#include <iostream>
#include <fstream>
#include "random.h"
#include "Metropolis.h"
#include "Metropolis2D.h"
#include "probDistr.h"
#include "probDensityTrial.h"
#include "BoltzmannWeight.h"
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

double Hpsi_over_psi(double, double, double);
double error(unsigned int, double, double);

int main(){

int distribution_type = 0; //0 means T is an uniform distribution, 1 means T is a gaussian distribution.
double mu = 0.80;
double sigma = 0.62;

probDensityTrial* myPDT = new probDensityTrial();

//Initial position: centre of the system.
Metropolis myMetropolis(myPDT, 0, distribution_type);

double position;
	
int nBlocks = 200;
int nInBlocks = 10000;
double mean_Energy;
double mean_Energy2;
double sum_Energy = 0;
double* new_muSigma;
int AccRej;

double old_energy, new_energy;

double AccAtt = 0.;

double sum = 0.;
double sum2 = 0.;

vector<double> estimations_Energy(nBlocks,0);
vector<double> estimations_uncertanties(nBlocks,0);

ofstream write_energy;
ofstream write_position;

write_energy.open("energy.dat");
write_position.open("positions.dat");

//set new parameters
myPDT->setParam(mu,sigma);

//Equilibration.
for(int i=0; i<1000; i++){
	myMetropolis.move();
}

//measure the new energy
for(int j=0; j<nBlocks; j++){

	sum_Energy = 0;
	
	for(int i=0; i<nInBlocks; i++){
		myMetropolis.move();
		position = myMetropolis.get_position();
		//cout << "position = " << position << endl;
		sum_Energy += Hpsi_over_psi(position, mu, sigma);
		write_position << position << endl;
	}
	
	//block results
	mean_Energy = sum_Energy/nInBlocks;
	mean_Energy2 = pow(mean_Energy,2);
	//cout << "mean Energy " << mean_Energy << endl;		

	sum += mean_Energy;
	sum2 += mean_Energy2;

	estimations_Energy[j] = sum/(j+1);

	//cout << "sum Energy = " <<estimations_Energy[j]<<endl;

	estimations_uncertanties[j] = error(j+1, sum/(j+1), sum2/(j+1));

	AccAtt = myMetropolis.get_Acc_rate();
	myMetropolis.reset_Acc_rate();
	//cout << "Acc/Att ratio = " << AccAtt << "   " << "Mean = " << estimations_Energy[j] << "  Error = " << estimations_uncertanties[j] << endl;	
	
	write_energy << setw(15) << estimations_Energy[j] << setw(15) << estimations_uncertanties[j] << endl;
	cout << setw(15) << estimations_Energy[j] << setw(15) << estimations_uncertanties[j] << endl;

}

write_energy.close();
write_position.close();

return 0;
}

double Hpsi_over_psi(double x, double mu, double sigma){

	double psi = exp(-pow(x-mu,2)/(2*sigma*sigma))+exp(-pow(x+mu,2)/(2*sigma*sigma));
	//cout << " psi = "<<exp(-pow(x-mu,2)/(2*sigma*sigma))<<endl;
	double Vpsi = (pow(x,4)-(5./2.)*pow(x,2))*psi;
	double Tpsi = -0.5*(1/(sigma*sigma))*((pow(x-mu,2)/(sigma*sigma)-1)*exp(-pow(x-mu,2)/(2*sigma*sigma))+(pow(x+mu,2)/(sigma*sigma)-1)*exp(-pow(x+mu,2)/(2*sigma*sigma))); //m, htagliato = 1
	/*if(contr == 0){
		cout << " x - mu ^ 2 = " << -pow(x-mu,2)<<endl;
		cout << " psi = "<<exp(-pow(x-mu,2)/(2*sigma*sigma))<<endl;
	}*/

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
