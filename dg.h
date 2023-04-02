#ifndef DG_H
#define DG_H 
#include "mesh.h"

namespace DG
{
	double calc_timestep(mesh &mesh1);
	void init_field(mesh &mesh1);
	void rhsboun(mesh &mesh1); // calculate the contribution of the domain integral to rhsel
	void rhsdomn(mesh &mesh1); //calculate the contribution of the boundary integral to rhsel
};	


//Flux methods 
namespace Flux
{
	
	// function to compute Roe's flux for a finite volume implementation.
	matrix2d fvRoe2d(mesh &mesh1,int i); //computes linearized flux using Roe's average pass mesh object by reference
	// here i is the index of the face we are calculating the flux for.
	matrix2d DGRoe2d(mesh &mesh1, int i);//flux for DG implementation 	
};


namespace ddt
{
	namespace explct
	{
		matrix2d RK3();
	};
};

#endif 

