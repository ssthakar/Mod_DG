#include "dgfem.h"
#include <chrono>
#include <iostream>
int main() 
{   
	std::ofstream file("out1.dat");
	mesh mesh1("channel.txt","control.txt");
	mesh1.preproc();
	std::cout<<mesh1.A.getval(1,1)<<std::endl;
	for(int i=0;i<mesh1.m_intface.size();i++)
	{
		for(int j=0;j<mesh1.m_intface[0].size();j++)
		{
			file<<mesh1.A.getval(i,j)<<" ";
		}
		file<<std::endl;
	}
    return 0;
}
