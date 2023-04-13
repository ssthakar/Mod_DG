#include "mesh.h"



//  unit test for the dg solver pre proc for now
int main()
{
	matrix2i A({1,1,1,1});
	printMatrix(A, "A.dat");
 	auto start = std::chrono::high_resolution_clock::now();
	grid::mesh mesh1("channel.txt","control.txt");
	grid::construct(mesh1);
	auto stop  = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop -start);
	std::cout<<duration.count()<<" microseconds "<<std::endl;
	return 0;
	
}
