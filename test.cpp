#include "mesh.h"
#include "vmatrix2.h"
#include "dg.h"


//  unit test for the dg solver pre proc for now
int main()
{
	auto start = std::chrono::high_resolution_clock::now();
	grid::mesh mesh1("mesh/channel.txt","mesh/control.txt");
	grid::construct(mesh1);
  printMatrix(mesh1.intface,"intface.out");
  printMatrix(mesh1.bface, "bface.out");
  printMatrix(mesh1.geoel,"geoel.out");
  printMatrix(mesh1.boun_geoface,"boun_geoface.out");
  printMatrix(mesh1.int_geoface, "int_geoface.out");
	soln soln1(mesh1,"mesh/control.txt");
  ddt::RK3::RK3_outer(mesh1,soln1);
  auto stop  = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop -start);
	std::cout<<duration.count()<<" microseconds "<<std::endl;
	feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
  return 0;
	
}
