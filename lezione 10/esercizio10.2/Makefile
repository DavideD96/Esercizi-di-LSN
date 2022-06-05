esercizio10_2: main.o random.o checkPath.o InitializePopulation.o cities.o geneticAlgorithm.o
	mpicxx main.o random.o checkPath.o InitializePopulation.o geneticAlgorithm.o cities.o -o esercizio10_2
main.o: main.cpp random.h checkPath.h InitializePopulation.h cities.h geneticAlgorithm.h
	mpicxx -c main.cpp -o main.o
random.o: random.cpp random.h
	mpicxx -c random.cpp -o random.o
checkPath.o: checkPath.cpp checkPath.h
	mpicxx -c checkPath.cpp -o checkPath.o
InitializePopulation.o: InitializePopulation.cpp InitializePopulation.h
	mpicxx -c InitializePopulation.cpp -o InitializePopulation.o
cities.o: cities.cpp cities.h
	mpicxx -c cities.cpp -o cities.o
geneticAlgorithm.o: geneticAlgorithm.cpp geneticAlgorithm.h
	mpicxx -c geneticAlgorithm.cpp -o geneticAlgorithm.o
