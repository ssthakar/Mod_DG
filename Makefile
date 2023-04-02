output: test.o mesh1.o 
	g++ mesh1.o test.o -o output

test.o: mesh1.cpp mesh.h
	g++ -c test.cpp

mesh1.o: mesh1.cpp mesh.h mMatrix.h mMatrix3.h 
	g++ -c mesh1.cpp

clean:
	rm *.o output
