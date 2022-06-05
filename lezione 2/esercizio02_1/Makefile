esercizio02_1: main.o random.o MCerror.o MCintegral.o
	g++ main.o random.o MCerror.o MCintegral.o -o esercizio02_1
main.o: main.cpp random.h
	g++ -c main.cpp -o main.o
random.o: random.cpp random.h
	g++ -c random.cpp -o random.o
MCerror.o: MCerror.cpp MCerror.h
	g++ -c MCerror.cpp -o MCerror.o
MCintegral.o: MCintegral.cpp MCintegral.h
	g++ -c MCintegral.cpp -o MCintegral.o
