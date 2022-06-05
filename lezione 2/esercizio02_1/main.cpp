#include <iostream>
#include <fstream>
#include "random.h"
#include "MCerror.h"
#include "MCintegral.h"
#include "cosine.h"
#include <cmath>
#include <iomanip>

using namespace std;

int main(){

//integrale con sampling uniforme
cosine* _cos_= new cosine; //funzione PI/2*cos(PI/2*x)
MCintegral myIntegral(_cos_, 0., 1., 100, 10000);
myIntegral.unifSampling();

double* risultatiMCUnif = myIntegral.getResults();
//for(int i=0; i<100; i++){
//	cout<<risultatiMC[i]<<endl;
//}
MCerror myMCerror1(risultatiMCUnif, 100);
myMCerror1.computeErrors();

double* medieUnif = myMCerror1.getMeans();
double* erroriUnif = myMCerror1.getErrors();

ofstream outData;
outData.open("unifSampl.dat");

for(int i=0; i<100; i++){
//	cout<<medieUnif[i]<<" "<<erroriUnif[i]<<endl;
	outData << setw(20) << medieUnif[i]<< setw(20) <<erroriUnif[i]<<endl;
}

outData.close();

//da completare...
//integrale con importance sampling usando una gaussiana centrata nell'origine (come scelgo la sigma??)
myIntegral.gaussSampling(0.47); //deviazione std
double* risultatiMCGauss = myIntegral.getResults();

MCerror myMCerror2(risultatiMCGauss, 100);
myMCerror2.computeErrors();

double* medieGauss = myMCerror2.getMeans();
double* erroriGauss = myMCerror2.getErrors();

outData.open("gaussSampl.dat");

for(int i=0; i<100; i++){
//	cout<<medieUnif[i]<<" "<<erroriUnif[i]<<"   "<<medieGauss[i]<<" "<<erroriGauss[i]<<endl;
	outData << setw(20) << medieGauss[i]<< setw(20) <<erroriGauss[i]<<endl;
}

outData.close();
/*
double singleRand;
double mean=0.;
double sigma=0.5;
double gaussVal;
double estimationGauss[100];
int i;

for(int j=0; j<100; j++){
	sumPart = 0.;
	i = 0;
	while(i<10000){
		singleRand = myRandom.Gauss(mean, sigma);
		if (singleRand>=0. && singleRand<1.){
			cout<<singleRand<<endl;
			sumPart += (M_PI/2.)*cos((M_PI/2.)*singleRand)/((1/sqrt(2*M_PI*pow(sigma,2)))*exp(-0.5*pow((singleRand-mean)/sigma,2)));
			i++;
		}
	}
	estimationGauss[j] = sumPart/10000;
	//cout<<estimationGauss[j]<<endl;
}
*/
return 0;
}
