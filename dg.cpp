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


//constructor for the flux object
FDS::RoeFlux::RoeFlux(matrix2d &Ul, matrix2d &Ur, double &nx, double &ny)
{

	assert(Ul.getncols() ==1 and Ur.getncols() ==1);
	double Pl = EOS::perf_gas(Ul);
	double Pr = EOS::perf_gas(Ur);
	avg.resize(5,1); //init container
	double Rij = sqrt(Ur(0,0)/Ul(0,0));
	avg(0,0) = Rij*Ul(0,0);
	avg(1,0) = (Rij*Ur(1,0)/Ur(0,0)+Ul(1,0)/Ul(0,0))/(1+Rij);
	avg(2,0) = (Rij*Ur(2,0)/Ur(0,0)+Ul(2,0)/Ul(0,0))/(1+Rij);
	avg(3,0) = (Rij*(Ur(3,0)+Pr)/Ur(0,0) + (Ul(3,0) + Pl)/Ul(0,0))/(Rij+1);
	avg(4,0) = (const_properties::gamma-1)*(avg(3,0) - 0.5*(avg(1,0)*avg(1,0) + avg(2,0)*avg(2,0)));
}

void FDS::RoeFlux::set_lamda()
{

}

void FDS::RoeFlux::set_W_amp()
{

}

void FDS::RoeFlux::flux_out()
{

}
