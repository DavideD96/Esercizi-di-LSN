esercizio05_1: main.o random.o Metropolis.o
	g++ main.o random.o Metropolis.o -o esercizio05_1
main.o: main.cpp Metropolis.h random.h probDistr.h groundState.h firstExcited.h
	g++ -c main.cpp -o main.o
random.o: random.cpp random.h
	g++ -c random.cpp -o random.o
Metropolis.o: Metropolis.cpp Metropolis.h
	g++ -c Metropolis.cpp -o Metropolis.o
