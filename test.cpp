#include "mesh.h"
#include "vmatrix2.h"
#include "dg.h"


//  unit test for the dg solver pre proc for now
int main()
{
  matrix3d L;
  std::ofstream file("unkel.out");
	auto start = std::chrono::high_resolution_clock::now();
	grid::mesh mesh1("mesh/naca.txt","mesh/control.txt");
	grid::construct(mesh1);
  std::cout<<"this the number of internal faces"<<mesh1.nintface<<std::endl;
  printMatrix(mesh1.intface,"intface.out");
  printMatrix(mesh1.bface, "bface.out");
  printMatrix(mesh1.boun_geoface,"boun_geoface.out");
  printMatrix(mesh1.esuel, "esuel.out");
	printMatrix(mesh1.int_geoface, "int_geoface.out");
  soln soln1(mesh1,"mesh/control.txt");
  //printMatrix(mesh1.U_infty,"initial.out");
  ddt::RK3::RK3_outer(mesh1,soln1);
  for(int i=0;i<mesh1.nelem;i++)
  {
     
    for(int j=0;j<mesh1.neqns;j++)
    {
      file<<mesh1.fvunkel(i,j)<<std::endl;//DG::U_at_poin(mesh1, mesh1.geoel(i, 2*j+5), mesh1.geoel(i, 2*j+6), i)(j,0)<<" ";
    }
    file<<std::endl;
  }
  file<<"wall y velocity"<<std::endl;
  file<<mesh1.geoel(3,6)<<std::endl;
  //42      0.14714821E+01      0.99689457E-01
//43      0.14968342E+01      0.99996180E-01
  double p1x = 0.14714821E+01;
  double p1y = 0.99689457E-01;
  double p2x = 0.14968342E+01 ;
  double p2y =  0.99996180E-01;
  double mag = len(p1x,p2x,p1y,p2y);
  double nx = (p2y-p1y)/mag;
  double ny = (p1x-p2x)/mag;
  matrix2d U = DG::U_at_poin(mesh1,mesh1.geoel(3,5),mesh1.geoel(3,6),3);
  file<<U(2,0)*ny  + U(1,0)*nx;
  printMatrix(mesh1.geoel,"geoel.out");
  
  auto stop  = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop -start);
	std::cout<<duration.count()<<" microseconds "<<std::endl;
	feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
  grid::post_proc::writevtk_mesh(mesh1,"this.vtk");
  return 0;
	
}
