#include <iostream>
#include <fstream>
#include <iomanip>
#include "random.h"
#include "Metropolis.h"
#include "probDistr.h"
#include "groundState.h"
#include "firstExcited.h"
#include <cmath>
#include <vector>

using namespace std;

//ATTENZIONE: AGGIUNGERE CONTROLLO DEL RATE DI ACCETTAZIONE

int main(){

ofstream Write100;

int distribution_type = 0; //0 means T is an uniform distribution, 1 means T is a gaussian distribution.

//Ground state | 100 >
groundState* myGS = new groundState;

//Initial position: centre of the atom (is surely the best choice).
Metropolis myMetropolis(myGS, 0,0,0);

//Equilibration.
double* position;

ofstream writeEq;
writeEq.open("radEq100origin.dat");

for(int i=0; i<1000; i++){
	myMetropolis.move(distribution_type); //0 means T is an uniform distribution, 1 means T is a gaussian distribution.
	position = myMetropolis.get_position();
	writeEq << i+1 << "  " << sqrt(position[0]*position[0]+position[1]*position[1]+position[2]*position[2]) << endl;
}

writeEq.close();

double sum_radius = 0;

int nBlocks = 200;
int nInBlocks = 500;
vector<double> mean_radius(nBlocks,0);
vector<double> mean_radius2(nBlocks,0);

Write100.open("100.txt");

for(int j=0; j<nBlocks; j++){

	sum_radius = 0;
	myMetropolis.reset_acc_tent();

	for(int i=0; i<nInBlocks; i++){
		myMetropolis.move(distribution_type); //0 means T is an uniform distribution, 1 means T is a gaussian distribution.
		position = myMetropolis.get_position();
		sum_radius += sqrt(position[0]*position[0]+position[1]*position[1]+position[2]*position[2]);
		//cout << " x = " <<position[0]<<"  "<<" y = " <<position[1]<<"  "<<" z = " <<position[2]<<endl;
		Write100 << position[0] << "  " << position[1] << "  " << position[2] << endl;
	}
	mean_radius[j] = sum_radius/nInBlocks;
	mean_radius2[j] = pow(mean_radius[j],2);

	cout << "Block " << j << ": acceptance rate" << " = " << myMetropolis.get_acc_rate() << endl;

}

Write100.close();

ofstream Radius100;
Radius100.open("radius100.txt");

vector<double> estimations_radius(nBlocks,0);
vector<double> estimations_uncertanties(nBlocks,0);

estimations_radius[0] = mean_radius[0];
double sum = mean_radius[0];
double sum2 = mean_radius2[0];

for(int i=1; i<nBlocks; i++){
	sum += mean_radius[i];
	sum2 += mean_radius2[i];
	
	estimations_radius[i] = sum/(i+1);
	estimations_uncertanties[i] = sqrt((1./((i+1)-1))*((1./(i+1))*sum2-pow(sum/(i+1),2)));

	//cout<<estimations_radius[i]<<"  "<<estimations_uncertanties[i]<<endl;
	Radius100 << i+1 << "  " << estimations_radius[i] << "  " << estimations_uncertanties[i] << endl;
	
}

Radius100.close();
cout << endl;
/*******************************************************************************************/

ofstream Write210;
Write210.open("210.txt");

//First excited | 210 >
firstExcited* myFE = new firstExcited;

//Initial position: centre of the atom.
Metropolis myMetropolis1(myFE, 0,0,-5);

writeEq.open("radEq210origin.dat");

//Equilibration.
for(int i=0; i<1000; i++){
	position = myMetropolis1.get_position();
	writeEq << i+1 << "  " << sqrt(position[0]*position[0]+position[1]*position[1]+position[2]*position[2]) << endl;
	myMetropolis1.move(distribution_type); //0 means T is an uniform distribution, 1 means T is a gaussian distribution.
}

writeEq.close();

double* position1;
double sum_radius1 = 0;

int nBlocks1 = 200;
int nInBlocks1 = 500;
vector<double> mean_radius1(nBlocks,0);
vector<double> mean_radius21(nBlocks,0);

for(int j=0; j<nBlocks1; j++){

	sum_radius1 = 0;
	myMetropolis1.reset_acc_tent();

	for(int i=0; i<nInBlocks1; i++){
		myMetropolis1.move(distribution_type); //0 means T is an uniform distribution, 1 means T is a gaussian distribution.
		position1 = myMetropolis1.get_position();
		sum_radius1 += sqrt(position1[0]*position1[0]+position1[1]*position1[1]+position1[2]*position1[2]);
		Write210 << position1[0] << "  " << position1[1] << "  " << position1[2] << endl;		
		//cout << " x = " <<position1[0]<<"  "<<" y = " <<position1[1]<<"  "<<" z = " <<position1[2]<<endl;
	}
	mean_radius1[j] = sum_radius1/nInBlocks1;
	mean_radius21[j] = pow(mean_radius1[j],2);

	cout << "Block " << j << ": acceptance rate" << " = " << myMetropolis1.get_acc_rate() << endl;

}

vector<double> estimations_radius1(nBlocks,0);
vector<double> estimations_uncertanties1(nBlocks,0);

estimations_radius1[0] = mean_radius1[0];
double sum1 = mean_radius1[0];
double sum21 = mean_radius21[0];

ofstream Radius210;
Radius210.open("radius210.txt");

for(int i=1; i<nBlocks1; i++){
	sum1 += mean_radius1[i];
	sum21 += mean_radius21[i];
	
	estimations_radius1[i] = sum1/(i+1);
	estimations_uncertanties1[i] = sqrt((1./((i+1)-1))*((1./(i+1))*sum21-pow(sum1/(i+1),2)));

	Radius210 <<i+1<<"  "<<estimations_radius1[i]<<"  "<<estimations_uncertanties1[i]<<endl;
	
}

Radius210.close();

return 0;
}
