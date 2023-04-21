#include "mesh.h"
#include "vmatrix2.h"
#include "dg.h"


//  unit test for the dg solver pre proc for now
int main()
{
	auto start = std::chrono::high_resolution_clock::now();
	grid::mesh mesh1("mesh/channel.txt","mesh/control.txt");
	grid::construct(mesh1);
  printMatrix(mesh1.bface, "bface.out");
	soln soln1(mesh1,"mesh/control.txt");
  ddt::RK3::RK3_outer(mesh1,soln1);
  printMatrix(mesh1.geoel, "geoel.out");
  auto stop  = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop -start);
	std::cout<<duration.count()<<" microseconds "<<std::endl;
	feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
  return 0;
	
}
