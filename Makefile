a.out: test.o mesh1.o dg.o mesh.h ddt.o
	g++   mesh1.o test.o dg.o ddt.o 

test.o: mesh1.cpp mesh.h dg.cpp dg.h ddt.cpp
	g++ -c test.cpp

mesh1.o: mesh1.cpp mesh.h vmatrix2.h mMatrix3.h ddt.cpp
	g++ -c mesh1.cpp 

dg.o: mesh1.cpp mesh.h dg.h ddt.cpp
	g++ -c dg.cpp

ddt.o: ddt.cpp dg.h mesh1.cpp mesh.h
	g++ -c ddt.cpp 

clean:
	rm *.o *.out 

