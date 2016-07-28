COMP=g++
FLAGS=-O3 
all:NOrMAL.o
	${COMP} ${FLAGS} NOrMAL.o main.cc -o NOrMAL 
NOrMAL.o: NOrMAL.h NOrMAL.cc
	${COMP} ${FLAGS} NOrMAL.cc -c 
test:all
	./NOrMAL test_config.txt DATA/Chr_01/T00f.txt DATA/Chr_01/T00r.txt test_results.txt	
