#include "mMatrix3.h"
#include "vmatrix2.h"
int main()
{
 
  
  vmatrix2<double> Ur(6,1);
  vmatrix2<double> Ul(4,1);
  print2Term(Ul);
  Ul = Ur;
  std::cout<<"ss \n";
  print2Term(Ul);
  mMatrix3<double> mat1;
  mat1.init(3, 4, 3);
  mat1.test_1(mat1);
  std::cout<<mat1(1,2,2);
  
}
