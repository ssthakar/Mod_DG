#include "mMatrix.h"
#include "mesh.h"
//  unit test for the dg solver pre proc for now
int main()
{	
	auto start = std::chrono::high_resolution_clock::now();
	grid::mesh mesh1("channel.txt","control.txt");
	grid::pre_proc::set_esup(mesh1);
	grid::pre_proc::set_esuel(mesh1);
	grid::pre_proc::set_intfafce(mesh1);
	grid::pre_proc::set_bface(mesh1);
	grid::pre_proc::set_geoel(mesh1);
	printMatrix(mesh1.esuel,"EEE.out");
	printMatrix(mesh1.intface,"fuckthisshit.out");
	printMatrix(mesh1.inpoel,"IIII");
	auto stop  = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop -start);
	std::cout<<duration.count()<<" microseconds "<<std::endl;
	std::cout<<mesh1.esuel(1371,0)<<mesh1.esuel(1371,1)<<mesh1.esuel(1371,2)<<std::endl;
	
	return 0;

}
