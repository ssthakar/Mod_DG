#include "dg.h"
#include "mesh.h"

//sub to push the contribution of the boundary integral to rhsel from bface only
void DG::rhsboun_bface(grid::mesh &mesh1)
{

}
//endsub


//sub to compute the contribution of the boundary integral at interfaces to RHSel from internal faces
void DG::rhsboun_iface(grid::mesh &mesh1)
{
	for(int i=0;i<mesh1.nintface;i++) //loop over all the internal faces 
	{
		for(int j=0;j<mesh1.ngauss_boun;j++) //loop over all gauss points
		{
			//get gauss coords
			//get conservative left variables at current gauss point 
			//get conservative right variables at current gauss poin 
			//store those conservative variables in a temp vector 
			// this vector is the function arg for Roes flux function
			//get flux at gauss coords from Roes Flux vector function
			//compute and push the contribution of the face to RHSel (add for left element and substract for right)
		}
	}
}
//endsub

// sub to push the contribution of domain integral to rhsel
void DG::rhsdomn(grid::mesh &mesh1)
{

}
//endsub 
