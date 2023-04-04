#ifndef DG_H
#define DG_H 
#include "mesh.h"

namespace DG
{
	void init_field(grid::mesh &mesh1);
	void rhsboun(grid::mesh &mesh1); // calculate the contribution of the domain integral to rhsel
	void rhsdomn(grid::mesh &mesh1); //calculate the contribution of the boundary integral to rhsel
};	


//Flux methods 
namespace Flux
{
	
	// function to compute Roe's flux for a finite volume implementation.
	matrix2d FVRoe(grid::mesh &mesh1,int i); //computes linearized flux using Roe's average pass mesh object by reference
	// here i is the index of the face we are calculating the flux for.
	matrix2d DGRoe(grid::mesh &mesh1, int i);//flux for DG implementation at face i 	
};


namespace ddt
{
	double calc_deltaT(grid::mesh &mesh1); //method to calculate delta T
	namespace explct
	{
		matrix2d fwd_euler(grid::mesh &mesh1, double &delta);
		matrix2d RK3(grid::mesh &mesh1,double &delta_t); // TVD runge Kutta 3rd order takes reference to mesh object and timestep 
	};
};

#endif 

