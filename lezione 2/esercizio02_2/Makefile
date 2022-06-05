esercizio02_2: main.o random.o cubicRW.o isotropicRW.o
	g++ main.o random.o cubicRW.o isotropicRW.o -o esercizio02_2
main.o: main.cpp random.h
	g++ -c main.cpp -o main.o
random.o: random.cpp random.h
	g++ -c random.cpp -o random.o
cubicRW.o: cubicRW.cpp cubicRW.h
	g++ -c cubicRW.cpp -o cubicRW.o
isotropicRW.o: isotropicRW.cpp isotropicRW.h
	g++ -c isotropicRW.cpp -o isotropicRW.o
clean:
	rm *.o
