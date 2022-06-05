esercizio01_1: main.o random.o distributionAnalyzer.o
	g++ main.o random.o distributionAnalyzer.o -o esercizio01_1
main.o: main.cpp random.h distributionAnalyzer.h
	g++ -c main.cpp -o main.o
random.o: random.cpp random.h
	g++ -c random.cpp -o random.o
distributionAnalyzer.o: distributionAnalyzer.cpp distributionAnalyzer.h
	g++ -c distributionAnalyzer.cpp -o distributionAnalyzer.o
