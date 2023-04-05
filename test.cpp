#include "mMatrix.h"
#include "mesh.h"
//  unit test for the dg solver pre proc for now
int main()
{	
	auto start = std::chrono::high_resolution_clock::now();
	grid::mesh mesh1("channel.txt","control.txt");
	grid::pre_proc::set_esup(mesh1);
	grid::pre_proc::set_esuel(mesh1);
	grid::pre_proc::set_bface(mesh1);
	grid::pre_proc::set_intface(mesh1);
	grid::pre_proc::set_geoel(mesh1);
	grid::pre_proc::set_massMat(mesh1);
	printMatrix(mesh1.geoel,"GGGG");
	auto stop  = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop -start);
	std::cout<<duration.count()<<" microseconds "<<std::endl;
	return 0;

}
