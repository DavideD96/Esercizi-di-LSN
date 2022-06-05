#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "priceComputer.h"

using namespace std;

int main(){

	double initialPrice = 100;
	double finalInstant = 1;
	double strikePrice = 100;
	double rf_interest_rate = 0.1;
	double volatility = 0.25;

	unsigned int nElements = 10000;
	unsigned int nBlocks = nElements/100;
	unsigned int elementInBlock = nElements/nBlocks;

	priceComputer myPricer(initialPrice,finalInstant,strikePrice,rf_interest_rate,volatility);

	//**************************************************
	//Prima parte: calcolo S(T).
	//**************************************************
	
	vector<double> callPrices(nElements,0);
	vector<double> putPrices(nElements,0);
	
	vector<double> meanCall(nBlocks,0);
	vector<double> meanPut(nBlocks,0);
	vector<double> sigmaCall(nBlocks,0);
	vector<double> sigmaPut(nBlocks,0);

	struct prices prices_;

	double meanC = 0;
	double meanP = 0;
	double meanC2 = 0;
	double meanP2 = 0;

	unsigned int index = 0;

	ofstream writeCall, writePut;
	writeCall.open("callPriceDir.dat");
	writePut.open("putPriceDir.dat");

	for(unsigned int j=0; j<nBlocks; j++){

		for(unsigned int i=0; i<elementInBlock; i++){
			
			prices_ = myPricer.callPutPrice();
			callPrices[index] = prices_.call;
			putPrices[index] = prices_.put;

			meanC = meanC + callPrices[index];
			meanP = meanP + putPrices[index];

			meanC2 = meanC2 + pow(callPrices[index],2);
			meanP2 = meanP2 + pow(putPrices[index],2);

			index++;
		}
		meanCall[j] = meanC/index;
		meanPut[j] = meanP/index;

		sigmaCall[j] = sqrt((meanC2/index-pow(meanC/index,2))/(index-1));
		sigmaPut[j] = sqrt((meanP2/index-pow(meanP/index,2))/(index-1));

		writeCall << setw(20) << meanCall[j] << setw(20) << sigmaCall[j] << endl;
		writePut << setw(20) << meanPut[j] << setw(20) << sigmaPut[j] << endl;
	}
	
	writeCall.close();
	writePut.close();

	//cout<<meanCall[999]<<endl;
	//cout<<sigmaCall[999]<<endl;

	//********************************************************
	//Seconda parte: calcolo S(T) discretizzando il cammino.
	//********************************************************

	vector<double> _callPrices(nElements,0);
	vector<double> _putPrices(nElements,0);
	
	vector<double> _meanCall(nBlocks,0);
	vector<double> _meanPut(nBlocks,0);
	vector<double> _sigmaCall(nBlocks,0);
	vector<double> _sigmaPut(nBlocks,0);

	struct prices _prices_;

	double _meanC = 0;
	double _meanP = 0;
	double _meanC2 = 0;
	double _meanP2 = 0;

	unsigned int _index = 0;

	unsigned int nSteps = 100;
	double deltaT = 1./nSteps;
	double S_;

	writeCall.open("callPriceDis.dat");
	writePut.open("putPriceDis.dat");

	for(unsigned int j=0; j<nBlocks; j++){

		for(unsigned int i=0; i<elementInBlock; i++){
			
			S_ = initialPrice;
			for(unsigned int k=0; k<nSteps; k++){
			
				S_ = myPricer.StStep(deltaT, S_);
			}

			double call_price = 0;
			double put_price = 0;
	
			if(S_-strikePrice > 0){
				call_price = exp(-rf_interest_rate*finalInstant)*(S_-strikePrice);
			}
			else{
				put_price = exp(-rf_interest_rate*finalInstant)*(strikePrice-S_);
			}
			_prices_ = myPricer.callPutPrice();
			callPrices[_index] = _prices_.call;
			putPrices[_index] = _prices_.put;

			_meanC = _meanC + callPrices[_index];
			_meanP = _meanP + putPrices[_index];

			_meanC2 = _meanC2 + pow(callPrices[_index],2);
			_meanP2 = _meanP2 + pow(putPrices[_index],2);

			_index++;
		}
		meanCall[j] = _meanC/_index;
		meanPut[j] = _meanP/_index;

		sigmaCall[j] = sqrt((_meanC2/_index-pow(_meanC/_index,2))/(_index-1));
		sigmaPut[j] = sqrt((_meanP2/_index-pow(_meanP/_index,2))/(_index-1));

		writeCall << setw(20) << meanCall[j] << setw(20) << sigmaCall[j] << endl;
		writePut << setw(20) << meanPut[j] << setw(20) << sigmaPut[j] << endl;		
	}	

writeCall.close();
writePut.close();

return 0;
}

