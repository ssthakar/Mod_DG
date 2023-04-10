#include "mMatrix.h"
#include "mesh.h"

class testclass
{
	public:
		int m;
		void compute(int nx);
};

void testclass::compute(int n)
{
	m = n*n;
}


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
	/*testclass shit;
	for(int i=0;i<10;i++)
	{
		shit.compute(i);
		std::cout<<shit.m<<std::endl;
	}
  */
	std::ofstream file("file_out");
	auto start = std::chrono::high_resolution_clock::now();
	grid::mesh mesh1("channel.txt","control.txt");
  grid::pre_proc::set_esup(mesh1);
  grid::pre_proc::set_geoel(mesh1);
  grid::pre_proc::set_esuel(mesh1);
	grid::pre_proc::set_bface(mesh1);
	grid::pre_proc::set_intface(mesh1);
	grid::pre_proc::set_geoel(mesh1);
	grid::pre_proc::set_massMat(mesh1);
	grid::pre_proc::set_int_geoface(mesh1);
	// write out stuff during cebuggin
  //printMatrix(mesh1.geoel,"GGGG");
	auto stop  = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop -start);
	std::cout<<duration.count()<<" microseconds "<<std::endl;

	return 0;
	
}
