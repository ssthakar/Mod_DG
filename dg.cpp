#include "dg.h" 
#include "mesh.h"



//TODO this is wrong implementation, there needs to be deltax and deltay in the denominator of the function
matrix2d DG::U_at_poin(grid::mesh &mesh1, double &gx, double &gy, int i) // function to get the conservative variable at any point withing an element
{
  matrix2d state;
  state.init(4, 1); // vector to return the state
  state(0, 0) = mesh1.unkel(i, 0, 0) + mesh1.unkel(i, 0, 1) * (gx - mesh1.geoel(i, 0)) + mesh1.unkel(i, 0, 2) * (gx - mesh1.geoel(i, 0));
  state(1, 0) = mesh1.unkel(i, 1, 0) + mesh1.unkel(i, 1, 1) * (gx - mesh1.geoel(i, 0)) + mesh1.unkel(i, 1, 2) * (gx - mesh1.geoel(i, 0));
  state(2, 0) = mesh1.unkel(i, 2, 0) + mesh1.unkel(i, 2, 1) * (gx - mesh1.geoel(i, 0)) + mesh1.unkel(i, 2, 2) * (gx - mesh1.geoel(i, 0));
  state(3, 0) = mesh1.unkel(i, 3, 0) + mesh1.unkel(i, 3, 1) * (gx - mesh1.geoel(i, 0)) + mesh1.unkel(i, 3, 2) * (gx - mesh1.geoel(i, 0));
  return state;
}


//function to get x-direction flux vector from state  
matrix2d DG::Fx(matrix2d &U) //used in domain integral only
{
	matrix2d fx(4,1);
	fx(0,0) = U(1,0); //rho U
  fx(1,0) = U(1,0)*U(1,0)/U(0,0)  + EOS::perf_gas(U); 
  fx(2,0) = U(1,0)*U(2,0)/U(0,0);
  fx(3,0) = (U(3,0)+EOS::perf_gas(U))*U(1,0)/U(0,0);
	return fx;
}
//function to get x-direction flux vector from state  
matrix2d DG::Fy(matrix2d &U) //used in domain integral only
{
	matrix2d fy(4,1);
	fy(0,0) = U(2,0); //rho U
  fy(1,0) = U(1,0)*U(2,0)/U(0,0);
  fy(2,0) = U(2,0)*U(2,0)/U(0,0)  + EOS::perf_gas(U); 
  fy(3,0) = (U(3,0)+EOS::perf_gas(U))*U(2,0)/U(0,0);
	return fy;
}


void DG::rhsdomn(grid::mesh &mesh1) //pass reference to the mesh object
{
  for(int i=0;i<mesh1.nelem;i++) //loop over all elems
  {
    for(int j=0;j<mesh1.ngauss_domn;j++) //loop over all gauss points
    {
      double gx,gy; //get gauss coorsd of current gaus point
      matrix2d U = DG::U_at_poin(mesh1, gx,gy,i); // get solution at current gauss point
      matrix2d fx,fy;
      fx = DG::Fx(U); //get flux function at current gauss point
      fy = DG::Fy(U);
      //get basis function derivatives at current gauss point
      double B2x,B2y,B3x,B3y;
      mesh1.rhsel(i,0,0) = fx(0,0)*B2x +fy(0,0)*B2y + mesh1.rhsel(i,0,0);
    }
  }
}


// sub to push the contribution of the boundary integral to rhsel from bface only
void DG::rhsboun_bface(grid::mesh &mesh1)
{
  assert(mesh1.unkel.size() > 1);        // make sure intialization of the flow field is done and the storage container is resized
  for (int i = 0; i < mesh1.nbface; i++) // loop over all boundary faces
  {
    // get the host cell number
    for (int j = 0; j < mesh1.ngauss_boun; j++)
    {
      // compute the numerical integral and push to rhsel location
    }
  }
}
// endsub


//TODO implement correct gauss coords rght now using nx and ny which is not correct 
// sub to compute the contribution of the boundary integral at interfaces to RHSel from internal faces
void DG::rhsboun_iface(grid::mesh &mesh1) 
{
	FDS::RoeFlux fluxobj; //instantiante flux object outside loop, to reduce the overhead of instantiating an object every time
  for (int i = 0; i < mesh1.nintface; i++) // loop over all the internal faces
  {
    double &nx = mesh1.int_geoface(i, 0); // xcomponent of area normal vector, needed for the flux
    double &ny = mesh1.int_geoface(i, 1); // y component of area normal vector
    int &le = mesh1.intface(i, 2);        // element to the left of the face when going from p1 to p2
    int &re = mesh1.intface(i, 3);        // element to the right of the face when going from p1 to p2
    for (int j = 0; j < mesh1.ngauss_boun; j++)
    {
			//get left and right states (conservative variable at the current gauss point)
			matrix2d Ul = DG::U_at_poin(mesh1,mesh1.int_geoface(i,j+2),mesh1.int_geoface(i,j+3),le);
			matrix2d Ur = DG::U_at_poin(mesh1,mesh1.int_geoface(i,j+2),mesh1.int_geoface(i,j+3),re);
      fluxobj.compute_flux(Ul, Ur, nx, ny); // compute the flux object using required data
      matrix2d flux = fluxobj.flux_intface();  // get the flux  at the interface
			mesh1.rhsel(le,0,0);
			mesh1.rhsel(le,0,1);
			mesh1.rhsel(le,0,2);
			//scatter the contributions to the rhsel data structure
		}
  }
}
// endsub


/*---------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//  imppmentation of Roes averaged flux, constructs flux object and has a method that returns the dissipative flux vector to be used in the rhsboun_iface subroutine
// instantiante object outside loop and then compute flux for every loop iteration
void FDS::RoeFlux::compute_flux(matrix2d &Ul, matrix2d &Ur, double &nx, double &ny) //compute function, computes the required data for other setter functions/kinda like a constructor
{
  assert(Ul.rows() == 4 and Ur.rows() == 4); // this implementation only works for 2d problems where the conservative vectors size is 4
  // compute and store in roe-averaged values //stores in other data that is needed for the caluclation of the dissipative flux
  Nx = nx;
  Ny = ny; // store in normal vectors for further usage
  assert(Ul.cols() == 1 and Ur.cols() == 1);
  double Pl = EOS::perf_gas(Ul); // compute left and right pressure states using EOS
  double Pr = EOS::perf_gas(Ur);
  deltaP = Pr - Pl; // compute delta values needed to compute wave amp
  deltaVx = Ur(1, 0) - Ul(1, 0);
  deltaVy = Ur(2, 0) - Ul(2, 0);
  delta_rho = Ur(0, 0) - Ul(0, 0);
  deltaVn = (Ur(1, 0) * Nx + Ur(2, 0) * ny) - (Ul(1, 0) * Nx + Ul(2, 0) * Ny);
  avg.init(5, 1); // init container , stores in roe agveraged vars
	avg_flux.init(4,1); //init container to store in central ave flux 
  double Rij = sqrt(Ur(0, 0) / Ul(0, 0));
  avg(0, 0) = Rij * Ul(0, 0);                                                                                            // rhoij
  avg(1, 0) = (Rij * Ur(1, 0) / Ur(0, 0) + Ul(1, 0) / Ul(0, 0)) / (1 + Rij);                                             // Vxij
  avg(2, 0) = (Rij * Ur(2, 0) / Ur(0, 0) + Ul(2, 0) / Ul(0, 0)) / (1 + Rij);                                             // Vyij
  avg(3, 0) = (Rij * (Ur(3, 0) + Pr) / Ur(0, 0) + (Ul(3, 0) + Pl) / Ul(0, 0)) / (Rij + 1);                               // Hij
  avg(4, 0) = sqrt((const_properties::gamma - 1) * (avg(3, 0) - 0.5 * (avg(1, 0) * avg(1, 0) + avg(2, 0) * avg(2, 0)))); // Cij speed of sound
	//compute averaged flux from right and left states at the cell interface
	avg_flux(0,0) = 0.5*(  Ul(1,0)*nx + Ul(2,0)* ny  + Ur(1,0)*nx + Ur(2,0)* ny);
	avg_flux(1,0) = 0.5*( (Ul(1,0)*Ul(1,0)+Pl)*nx + Ul(1,0)*Ul(2,0)/Ul(0,0)*ny + (Ur(1,0)*Ur(1,0)+Pr)*nx + Ur(1,0)*Ur(2,0)/Ur(0,0)*ny);
	avg_flux(2,0) = 0.5*(Ul(1,0)*Ul(2,0)/Ul(0,0)*nx + (Ul(2,0)*Ul(2,0)/Ul(0,0)+Pl)*ny + Ur(1,0)*Ur(2,0)/Ur(0,0)*nx + (Ur(2,0)*Ur(2,0)/Ur(0,0)+Pr)*ny);
	avg_flux(3,0) = 0.5*((Ul(3,0)+Pl)*Ul(1,0)/Ul(0,0)*nx + (Ul(3,0)+Pl)*Ul(2,0)/Ul(0,0)*ny + (Ur(3,0)+Pr)*Ur(1,0)/Ur(0,0)*ny+ (Ur(3,0)+Pr)*Ur(2,0)/Ur(0,0)*ny);
}

// compute the eigen values of the A matrix
void FDS::RoeFlux::set_lamda()
{
  // first two lamdas are the same
  lamda(3, 1); // init lamda
  lamda(0, 0) = avg(1, 0) * Nx + avg(2, 0) * Ny;
  lamda(1, 0) = lamda(0, 0) + avg(4, 0);
  lamda(2, 0) = lamda(0, 0) - avg(4, 0);
}

// compute the wave amplitudes and store them
void FDS::RoeFlux::set_W_amp()
{
  W_amp(4, 0); // inti container and compute store in wave amplitudes
  W_amp(0, 0) = delta_rho - deltaP / pow(avg(4, 0), 2);
  W_amp(1, 0) = deltaVx * Ny - deltaVy * Nx;
  W_amp(2, 0) = deltaVn + deltaP / (avg(0, 0) * avg(4, 0));
  W_amp(3, 0) = -1 * deltaVn + deltaP / (avg(0, 0) * avg(4, 0));
}
//compute the dissipative flux that is substracted from the central averaged flux at a gauss point
void FDS::RoeFlux::set_diss_flux()
{
  double r = avg(0, 0) / (2 * avg(4, 0));

  dissFlux.init(4, 1);
  // this is the dissipative flux that is substracted from the central average
  double a0 = std::fabs(lamda(0, 0)) * W_amp(0, 0);
  double a1 = std::fabs(lamda(0, 0)) * W_amp(0, 0);
  double a2 = std::fabs(lamda(1, 0)) * W_amp(1, 0);
  double a3 = std::fabs(lamda(2, 0)) * W_amp(2, 0);
  dissFlux(0, 0) = (a0)*1 + 0.5 * avg(0, 0) / avg(4, 0) * lamda(2, 0) * +0.5 * lamda(3, 0) * avg(0, 0) / avg(4, 0);
  dissFlux(1, 0) = a0 * avg(1, 0) + a1 * avg(4, 0) * Ny + a2 * r * (avg(1, 0) + avg(4, 0) * Nx) + a3 * r * (avg(1, 0) - avg(4, 0) * Nx);
  dissFlux(2, 0) = a0 * avg(2, 0) - a1 * avg(4, 0) * Nx + a2 * r * (avg(2, 0) + avg(4, 0) * Nx) + a3 * r * (avg(2, 0) - avg(4, 0) * Nx);
  dissFlux(3, 0) = a0 * 0.5 * (avg(1, 0) * avg(1, 0) + avg(2, 0) * avg(2, 0)) + a1 * avg(4, 0) * (avg(1, 0) * Ny - avg(2, 0) * Nx) + a2 * (avg(3, 0) + avg(4, 0) * (avg(1, 0) * Nx + avg(2, 0) * Ny)) + a3 * (avg(3, 0) - avg(4, 0) * (avg(1, 0) * Nx + avg(2, 0) * Ny));
}
 //get interface flux by calling this function in the rhsboun subs
matrix2d FDS::RoeFlux::flux_intface()
{
	matrix2d intface_flux(4,1);
	set_lamda();
  set_W_amp();
  set_diss_flux();
	intface_flux(0,0)  = avg_flux(0,0) - dissFlux(0,0);
	intface_flux(1,0)  = avg_flux(1,0) - dissFlux(1,0);
	intface_flux(2,0)  = avg_flux(2,0) - dissFlux(2,0);
	intface_flux(3,0)  = avg_flux(3,0) - dissFlux(3,0);
  return intface_flux; // returns the dissipative flux , this is not the flux function
}
// end class definition



//function to calculate local time step for cell i 
double ddt::local_ts(grid::mesh &mesh1, int &i)
{
	double delta_T;
	

	return delta_T;
}


