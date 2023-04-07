#include "mMatrix.h"
#include "mesh.h"
std::vector<double> testFunc()
{
	std::vector<double> B;
	B.resize(3);
	B[0] = 1;
	B[1] = B[0]+6;
	B[2] = B[1]+B[0]*140;
	return B;
}
//  unit test for the dg solver pre proc for now
int main()
{
	std::vector<double> Temp = testFunc();
	std::ofstream file("file_out");
	auto start = std::chrono::high_resolution_clock::now();
	grid::mesh mesh1("channel.txt","control.txt");
	grid::pre_proc::set_esup(mesh1);
	grid::pre_proc::set_esuel(mesh1);
	grid::pre_proc::set_bface(mesh1);
	grid::pre_proc::set_intface(mesh1);
	grid::pre_proc::set_geoel(mesh1);
	grid::pre_proc::set_massMat(mesh1);
	grid::pre_proc::set_int_geoface(mesh1);
	// write out stuff durint debugging
	/*for(int i=0;i<mesh1.nelem;i++)
	{
		file<<(mesh1.geoel(i,1)*mesh1.geoel(i,3) - mesh1.geoel(i,2)*mesh1.geoel(i,2))<<std::endl;
	} */
	printMatrix(mesh1.geoface,"GGGG");
	auto stop  = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop -start);
	std::cout<<duration.count()<<" microseconds "<<std::endl;
	for(int i=0;i<Temp.size();i++)
	{
		std::cout<<Temp[i]<<std::endl;
	}

	return 0;
	
}
