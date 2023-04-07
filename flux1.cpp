#include "dg.h"
#include "mesh.h"

matrix2d Flux::Roe::intface_flux(grid::mesh &mesh1, int i)
{
	matrix2d flux(4,1);
	int &L = mesh1.intface(i,2);
	int &R = mesh1.intface(i,3);
	matrix2d UL = 
	return flux;
}
