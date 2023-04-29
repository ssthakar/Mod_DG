#include "mMatrix3.h"

int main()
{
  mMatrix3<double> mat1;
  mat1.init(3, 4, 3);
  mat1.test_1(mat1);
  std::cout<<mat1(1,2,2);
  
}
