#include "mesh.h"
#include "dg.h"


std::vector<double> Flux::Roe::Cons2Prim(std::vector<double> &state)
{
	std::vector<double> Prim;
	Prim.resize(5); // add one place holder for pressure computation 
	Prim[0] = state[0]; //density
	Prim[1] = state[1]/state[0]; //u or x vel 
	Prim[2] = state[2]/state[0]; //v or y vel 
	Prim[3] = state[3]/state[0]; // total energy e
	Prim[4] = EOS::perf_gas(state); // calculate and store pressure 
	return Prim;
}


//function to compute the Roe_averaged variables 
std::vector<double> Flux::Roe::Roe_avg(std::vector<double> &L_prim, std::vector<double> &R_prim) //take in prim vars 
{
	std::vector<double> Ravg;
	Ravg.resize(4); // this is only for 2d not considering general case 
	double hL = (L_prim[0]*L_prim[3] + L_prim[4])/L_prim[0]; //calculate the enthalpy for the left side
	double hr = (R_prim[0]*R_prim[3] + R_prim[4])/R_prim[0]; //calculate the enthalpy for the right side 
	return Ravg;
}


//function to compute the wave speed 
std::vector<double> Flux::Roe::wave_speed(std::vector<double> &vec)
{
	std::vector<double> ws;
	return ws;
}

//function to compute the wave strength given 
std::vector<double> Flux::Roe::wave_strength(std::vector<double> &vec)
{
	std::vector<double> wstrngth;
	return wstrngth;
}
// function to compute the Roes average flux for ith internal face
// takes in comnined vector of left and right conservative variables to compute the flux at Guass point
//                       rhoL rhoUL rhoVL rhoEL rhoR  rhoUR	 rhoVR 	rhoER
//format for arg vector (vl1  vl2   vl3   vl4   vr1   vr2    vr3    vr4) //conservative variable vector
//											 0		1	 	  2		  3	    4		  5		   6		  7	
//
matrix2d Flux::Roe::DGRoe2d(grid::mesh &mesh1,std::vector<double> &LR)
{
	//assert(args.size() == 2*mesh1.neqns); //make sure the right vector is passed toa the function
	matrix2d flux(mesh1.neqns,1); //init empty container to store in flux at interface
	//compute roe averaged quantities from conservative vars
	std::vector<double> Left_state(LR.begin(),LR.begin()+4);
	std::vector<double> Right_state(LR.begin()+4, LR.end());
	std::vector<double> PrimL = Flux::Roe::Cons2Prim(Left_state); //get primary variables on the left  
	std::vector<double> PrimR = Flux::Roe::Cons2Prim(Right_state); // get primary variables on the right
	//compute enthalpy for use in Roes averages
	//compute Roes average vars
	//compute wave speeds
	//compute wave strengths
	// for loop to compute flux at interface 
	//
	
	return flux;

}




