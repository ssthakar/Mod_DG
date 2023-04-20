#include "mesh.h"
#include "vmatrix2.h"
#include "dg.h"


//  unit test for the dg solver pre proc for now
int main()
{
	matrix2d A(4,1);
	std::cout<<A.rows()<<std::endl;
	auto start = std::chrono::high_resolution_clock::now();
	grid::mesh mesh1("mesh/test.txt","mesh/control.txt");
	grid::construct(mesh1);
	printMatrix(mesh1.esuel, "esuel.dat");
	soln soln1(mesh1,"mesh/control.txt");
	std::cout<<mesh1.unkel(3,2,0)<<std::endl;
  ddt::RK3::RK3_outer(mesh1,soln1);
	printMatrix(mesh1.geoel,"geoel.out");
	auto stop  = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop -start);
	std::cout<<duration.count()<<" microseconds "<<std::endl;
	std::cout<<mesh1.unkel(3,2,0)<<std::endl;
	return 0;
	
}
