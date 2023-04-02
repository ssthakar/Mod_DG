#include "mMatrix.h"
#include "mesh.h"
//  unit test for the dg solver pre proc for now
int main()
{	
	auto start = std::chrono::high_resolution_clock::now();
	mesh mesh1("channel.txt","control.txt");
	mesh1.pre_proc();
	mMatrix2<int> test(mesh1.inpoel);
	matrix3d test1;
	int & shri = mesh1.intface(1,1);
	test1.resize(5,5,5);
	std::cout<<test1(1,1,1)<<" this shit better work"<<std::endl;
	printMatrix(mesh1.bface,"bface.out");
	printMatrix(mesh1.intface,"intface.out");
	auto stop  = std::chrono::high_resolution_clock::now();
	printMatrix(mesh1.coords,"out2.dat");//std::cout<<mesh1.inpoel.getval(6,1);
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop -start);
	std::cout<<duration.count()<<" microseconds "<<std::endl;
	//std::cout<<mesh1.npoin<<" "<<mesh1.esup2(mesh1.npoin,0)<<std::endl;
	std::cout<<shri;
	return 0;

}
