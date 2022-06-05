#ifndef _distributionAnalyzer_h_
#define _distributionAnalizer_h_

class distributionAnalyzer{

public:
	distributionAnalyzer(double*, unsigned int); //costruttore che inizializza _data con l'array di dati estratti (ad esempio n° tra 0 e 1). (dati, dimensione array)
	void setExtrema(double, double); //assegna gli estremi della distribuzione (estremo sx, estremo dx).
	void setPeriod(unsigned int);	//assegna il numero di dati in un blocco.
	void setNumBlock(unsigned int); //assegna il numero di blocchi.
	void setIntervals(unsigned int); //assegna il nuumero di sottointervalli in cui dividere il dominio.
	int** getCountsInInt();
	void dataInIntervals(); //calcola i dati in sottointervalli uguali (usa la funzione set per stabilire il numero di intervalli e il numero di blocchi!).
	double* chi2(); //calcola il chi2 per ogni blocco di numeri, resituisce un array con tutti i chi2 calcolati.
	~distributionAnalyzer();
	
private:
	double* _randNumbers;
	unsigned int _M; //numero dati
	double _extSx;
	double _extDx;
	unsigned int _K; //numero intervallini //INUTILI?
	unsigned int _period; //numero dati in un blocco
	unsigned int _N; //numero blocchi (dipendente da _period e _M)
	int** _counts_per_bin;
	
};

#endif
