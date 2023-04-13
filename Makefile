a.out: test.o mesh1.o dg.o mesh.h
	g++ -O1  mesh1.o test.o 

test.o: mesh1.cpp mesh.h
	g++ -c test.cpp

mesh1.o: mesh1.cpp mesh.h mMatrix.h mMatrix3.h 
	g++ -c mesh1.cpp mesh.h

dg.o: mesh1.cpp mesh.h mMatrix3.h mMatrix.h dg.h
	g++ -c dg.cpp

clean:
	rm *.o *.dat a.out
