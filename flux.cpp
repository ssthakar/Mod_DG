#include "mesh.h"

//format for intface row for a face with index i 
// (ip1 | ip2 |cell on left | cell on right)
matrix2d Flux::fvRoe2d(mesh &mesh1, int i)
{
	matrix2d flux(mesh1.neqns);//init matrix to store in flux, return this 
	assert(mesh1.ndegr == 0); //this function is only for finite volume 
	//get the element numbers on the right and left intface 
	int &L = mesh1.intface(i,2);
	int &R  = mesh1.intface(i,3);
	//format of solution storage for element L or R 
	// fv_U(i,) = rho | U | V | E 
	// declare variables needed to compute flux
	//double rho_R,u_L,u_R,v_L,v_R,E_L,E_R,p;
	double & rho_L = mesh1.fv_U(L,0); //density left state
	double & rho_R = mesh1.fv_U(R,0);
	double & u_L = mesh1.fv_U(L,1); //x component left state 
	double & u_R = mesh1.fv_U(R,1);
	double & v_L = mesh1.fv_U(L,2);
	double & v_R = mesh1.fv_U(R,2);
	double & E_L = mesh1.fv_U(L,3); //energy left state 
	double & E_R = mesh1.fv_U(R,3);
	//compute the pressure for both states of the solution;
	double p_L = EOS::perf_gas(rho_L,E_L,u_L,v_L);
	double p_R = EOS::perf_gas(rho_R,E_R,u_R,v_R);
	double R_LR = sqrt(rho_L/rho_R);
	double rho_LR = R_LR*rho_L;
	
	return flux;

}




