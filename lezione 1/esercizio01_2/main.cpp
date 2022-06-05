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

//Write on file
ofstream WriteUnif;
ofstream WriteExp;
ofstream WriteLor;

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

double* expNumbers = new double[M]; 
double* LorentzNumbers = new double[M]; 
double* unifNumbers = new double[M];

WriteUnif.open("unifS_1.txt");
WriteExp.open("expS_1.txt");
WriteLor.open("lorS_1.txt");

//carico vettore con numeri estratti dalla distribuzione uniforme.
//carico vettore con numeri estratti dalla distribuzione esponenziale con lambda = 1.
//carico vettore con numeri estratti dalla distribuzione lorentziana con gamma = 1, media = 0;

for (int i=0; i<M; i++){
	unifNumbers[i] = myRandom.Rannyu(0., 1.);
	expNumbers[i] = myRandom.Exp(1.);
	LorentzNumbers[i] = myRandom.Lorentz(0., 1.);

	WriteUnif<<unifNumbers[i]<<endl;
	WriteExp<<expNumbers[i]<<endl;
	WriteLor<<LorentzNumbers[i]<<endl;
	
}

WriteUnif.close();
WriteExp.close();
WriteLor.close();

//S_2

WriteUnif.open("unifS_2.txt");
WriteExp.open("expS_2.txt");
WriteLor.open("lorS_2.txt");

unsigned int Ns = 2;
unsigned int N = M/Ns;

double* S_2unif = new double[N];
double* S_2exp = new double[N];
double* S_2Lor = new double[N];

for (int i=0; i<N; i++){
	for (int j=0; j<Ns; j++){
		S_2unif[i] += unifNumbers[Ns*i+j];
		S_2exp[i] += expNumbers[Ns*i+j];
		S_2Lor[i] += LorentzNumbers[Ns*i+j];
	}
	WriteUnif<<S_2unif[i]<<endl;
	WriteExp<<S_2exp[i]<<endl;
	WriteLor<<S_2Lor[i]<<endl;
}

WriteUnif.close();
WriteExp.close();
WriteLor.close();

//S_10

WriteUnif.open("unifS_10.txt");
WriteExp.open("expS_10.txt");
WriteLor.open("lorS_10.txt");

Ns = 10;
N = M/Ns;

double* S_10unif = new double[N];
double* S_10exp = new double[N];
double* S_10Lor = new double[N];

for (int i=0; i<N; i++){
	for (int j=0; j<Ns; j++){
		S_10unif[i] += unifNumbers[Ns*i+j];
		S_10exp[i] += expNumbers[Ns*i+j];
		S_10Lor[i] += LorentzNumbers[Ns*i+j];
	}
	WriteUnif<<S_10unif[i]<<endl;
	WriteExp<<S_10exp[i]<<endl;
	WriteLor<<S_10Lor[i]<<endl;
}

WriteUnif.close();
WriteExp.close();
WriteLor.close();

//S_100

WriteUnif.open("unifS_100.txt");
WriteExp.open("expS_100.txt");
WriteLor.open("lorS_100.txt");

Ns = 100;
N = M/Ns;

double* S_100unif = new double[N];
double* S_100exp = new double[N];
double* S_100Lor = new double[N];

for (int i=0; i<N; i++){
	for (int j=0; j<Ns; j++){
		S_100unif[i] += unifNumbers[Ns*i+j];
		S_100exp[i] += expNumbers[Ns*i+j];
		S_100Lor[i] += LorentzNumbers[Ns*i+j];
	}
	WriteUnif<<S_100unif[i]<<endl;
	WriteExp<<S_100exp[i]<<endl;
	WriteLor<<S_100Lor[i]<<endl;
}

WriteUnif.close();
WriteExp.close();
WriteLor.close();

}
