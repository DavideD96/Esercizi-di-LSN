/*
INTESTAZIONE
*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "random.h"

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

unsigned int M = 1000000; //un milione di tiri
unsigned int nBlk = 100;
unsigned int nPerBlock = M/nBlk;

double L = 1.; //lunghezza bastoncino
double d = 2.;
double pos1, pos2;

double cross[nBlk][2];
cross[0][1] = 0; //incertezza primo blocco
int cont;
double avg = 0;
double avg2 = 0;

ofstream outBuff;
outBuff.open("Buffon.dat");

for(int i=0; i<nBlk; i++){

	cont = 0;
	for(int j=0; j<nPerBlock; j++){
		pos1 = d*myRandom.Rannyu();
		pos2 = pos1 + L*sin(myRandom.Rannyu(0,M_PI/2));
		if(pos2 > d){
			cont++;
		}
	}
	avg += 2.*L/(((double)cont/nPerBlock)*d);
	avg2 += pow(2.*L/(((double)cont/nPerBlock)*d),2);
	cross[i][0] = avg/(i+1);
	if(i != 0){
		cross[i][1] = sqrt((1./i)*(avg2/(i+1)-cross[i][0]*cross[i][0]));
		//cout << avg << "  " << avg2 << endl;
	}
	outBuff << cross[i][0] << " " << cross[i][1] << endl;
}

outBuff.close();

return 0;
}
