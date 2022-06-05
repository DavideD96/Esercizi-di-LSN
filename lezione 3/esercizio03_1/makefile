esercizio03_1: main.o random.o priceComputer.o
	g++ main.o random.o priceComputer.o -o esercizio03_1
main.o: main.cpp random.h
	g++ -c main.cpp -o main.o
random.o: random.cpp random.h
	g++ -c random.cpp -o random.o
priceComputer.o: priceComputer.cpp priceComputer.h
	g++ -c priceComputer.cpp -o priceComputer.o
