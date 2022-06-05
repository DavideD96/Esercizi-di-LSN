esercizio10_1: main.o random.o InitializePopulation.o cities.o SimulatedAnnealing.o
	g++ main.o random.o InitializePopulation.o SimulatedAnnealing.o cities.o -o esercizio10_1
main.o: main.cpp random.h InitializePopulation.h cities.h SimulatedAnnealing.h
	g++ -c main.cpp -o main.o
random.o: random.cpp random.h
	g++ -c random.cpp -o random.o
InitializePopulation.o: InitializePopulation.cpp InitializePopulation.h
	g++ -c InitializePopulation.cpp -o InitializePopulation.o
cities.o: cities.cpp cities.h
	g++ -c cities.cpp -o cities.o
SimulatedAnnealing.o: SimulatedAnnealing.cpp SimulatedAnnealing.h
	g++ -c SimulatedAnnealing.cpp -o SimulatedAnnealing.o
