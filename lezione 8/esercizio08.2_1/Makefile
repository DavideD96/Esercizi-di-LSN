esercizio08_2: main.o random.o Metropolis.o Metropolis2D.o
	g++ main.o random.o Metropolis.o Metropolis2D.o -o esercizio08_2
main.o: main.cpp Metropolis.h Metropolis2D.h random.h probDistr.h probDensityTrial.h BoltzmannWeight.h
	g++ -c main.cpp -o main.o
random.o: random.cpp random.h
	g++ -c random.cpp -o random.o
Metropolis.o: Metropolis.cpp Metropolis.h
	g++ -c Metropolis.cpp -o Metropolis.o
Metropolis2D.o: Metropolis2D.cpp Metropolis2D.h
	g++ -c Metropolis2D.cpp -o Metropolis2D.o
