#include "priceComputer.h"

priceComputer::priceComputer(double ip, double fi, double sp, double ir, double vol){

	_initialPrice = ip;
	_finalInstant = fi;
	_strikePrice = sp;
	_rf_interest_rate = ir;
	_volatility = vol;

	ifstream inData;

	int seed[4];
	int p1, p2;

	inData.open("Primes.txt");

	if (inData.is_open()){
		inData >> p1 >> p2;
	} 
	else{
		std::cerr << "Unable to open Primes" <<std::endl;
	}
	inData.close(); 

	inData.open("seed.in");
	string property;

	if (inData.is_open()){
		while(!inData.eof()){
			inData >> property;
			if(property == "RANDOMSEED"){
				inData >> seed[0] >> seed[1] >> seed[2] >>seed[3];
			_gen.SetRandom(seed, p1, p2);
			}
		}
		inData.close();
	}
}

struct prices priceComputer::callPutPrice(){

	double randValue = _gen.Gauss(0,_finalInstant);
	double St = _initialPrice*exp((_rf_interest_rate-0.5*pow(_volatility,2))*_finalInstant+_volatility*randValue);
	double call_price = 0;
	double put_price = 0;	
	if(St-_strikePrice > 0){
		call_price = exp(-_rf_interest_rate*_finalInstant)*(St-_strikePrice);
	}
	else{
		put_price = exp(-_rf_interest_rate*_finalInstant)*(_strikePrice-St);
	}
	/*double d1 = (1./(_volatility*sqrt(_finalInstant)))*(log(St/_strikePrice)+(_rf_interest_rate+pow(_volatility,2)*_finalInstant/2.));
	double d2 = d1 -_volatility*sqrt(_finalInstant);
	double Nd1 = 0.5*(1+erf(d1/(sqrt(2.))));
	double Nd2 = 0.5*(1+erf(d2/(sqrt(2.))));

	double price = St*Nd1-_strikePrice*exp(_rf_interest_rate*_finalInstant)*Nd2;*/
	struct prices pr;
	pr.call = call_price;
	pr.put = put_price;

	return pr;
}

double priceComputer::StStep(double deltaT, double previousSt){	

	double randValue = _gen.Gauss(0,1);
	double St = previousSt*exp((_rf_interest_rate-0.5*pow(_volatility,2))*deltaT+_volatility*randValue*sqrt(deltaT));

	return St;
}
