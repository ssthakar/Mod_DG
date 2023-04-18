#include "mesh.h"
#include "vmatrix2.h"
#include "dg.h"


//  unit test for the dg solver pre proc for now
int main()
{

	auto start = std::chrono::high_resolution_clock::now();
	grid::mesh mesh1("mesh/naca.txt","mesh/control.txt");
	grid::construct(mesh1);
	soln soln1(mesh1,"mesh/control.txt");
	
	auto stop  = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop -start);
	std::cout<<duration.count()<<" microseconds "<<std::endl;
	return 0;
	
}
