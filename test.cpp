#include "mesh.h"
#include "vmatrix2.h"



//  unit test for the dg solver pre proc for now
int main()
{

 	auto start = std::chrono::high_resolution_clock::now();
	grid::mesh mesh1("mesh/test.txt","mesh//control.txt");
	grid::construct(mesh1);
  printMatrix(mesh1.geoel, "out/testout.dat");
	auto stop  = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop -start);
	std::cout<<duration.count()<<" microseconds "<<std::endl;
	return 0;
	
}
