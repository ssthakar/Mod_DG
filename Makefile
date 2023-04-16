a.out: test.o mesh1.o dg.o mesh.h ddt.o
	g++ -O1  mesh1.o test.o 

test.o: mesh1.cpp mesh.h
	g++ -c test.cpp

mesh1.o: mesh1.cpp mesh.h vmatrix2.h mMatrix3.h 
	g++ -c mesh1.cpp mesh.h

dg.o: mesh1.cpp mesh.h mMatrix3.h vmatrix2.h dg.h ddt.cpp
	g++ -c dg.cpp

ddt.o: ddt.cpp dg.h mesh1.cpp mesh.h vmatrix2.h mMatrix3.h
	g++ -c ddt.cpp dg.h

clean:
	rm *.o  a.out

