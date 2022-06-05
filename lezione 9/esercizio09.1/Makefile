esercizio09_1: main.o random.o checkPath.o InitializePopulation.o cities.o geneticAlgorithm.o
	g++ main.o random.o checkPath.o InitializePopulation.o geneticAlgorithm.o cities.o -o esercizio09_1
main.o: main.cpp random.h checkPath.h InitializePopulation.h cities.h geneticAlgorithm.h
	g++ -c main.cpp -o main.o
random.o: random.cpp random.h
	g++ -c random.cpp -o random.o
checkPath.o: checkPath.cpp checkPath.h
	g++ -c checkPath.cpp -o checkPath.o
InitializePopulation.o: InitializePopulation.cpp InitializePopulation.h
	g++ -c InitializePopulation.cpp -o InitializePopulation.o
cities.o: cities.cpp cities.h
	g++ -c cities.cpp -o cities.o
geneticAlgorithm.o: geneticAlgorithm.cpp geneticAlgorithm.h
	g++ -c geneticAlgorithm.cpp -o geneticAlgorithm.o
