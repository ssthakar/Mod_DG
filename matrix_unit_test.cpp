#include "mMatrix.h"

int main()
{
	mMatrix2<double> A(3,1);
	A(0,0) = 1;
	A(1,0) = 2;
	A(2,0) = 3;
	mMatrix2<double> B = A*A;
	printMatrix(B, "matrix.out");
}
