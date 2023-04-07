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


//function to compute the Roe_averaged variables taken in vector of primitive vars and pressure
std::vector<double> Flux::Roe::Roe_avg(std::vector<double> &L_prim, std::vector<double> &R_prim) //take in prim vars 
{
	std::vector<double> Ravg;
	Ravg.resize(5); // this is only for 2d not considering general case 
	double hL = (L_prim[0]*L_prim[3] + L_prim[4])/L_prim[0]; //calculate the enthalpy for the left side
	double hR = (R_prim[0]*R_prim[3] + R_prim[4])/R_prim[0];//calculate the enthalpy for the right side 
	double Rij = sqrt(R_prim[0]/L_prim[0]); //sqrt(rhoR/rhoL)
	Ravg[0] = Rij*L_prim[0]; //rho_ij
	Ravg[1] = (Rij*R_prim[1] + L_prim[1])/(Rij+1); //uij
	Ravg[2] = (Rij*R_prim[2] + L_prim[2])/(Rij+1); //vij
	Ravg[3] = (Rij*hR+hL)/(Rij+1); //Hij
	Ravg[4] = sqrt((const_properties::gamma-1)*(Ravg[3] - 0.5*(Ravg[1]*Ravg[1]+Ravg[2]*Ravg[2]))); //cij 
	return Ravg;
}


//function to compute the wave speed takes in vector of roe averaged variables 
std::vector<double> Flux::Roe::wave_speed(std::vector<double> &ravg_vector, double &nx, double &ny)
{
	std::vector<double> ws;
	ws.resize(3);
	ws[0] =ravg_vector[1]*nx+ravg_vector[2]*ny;
	ws[1] = ws[0]+ravg_vector[4];
	ws[2] = ws[0] + ravg_vector[4];
	return ws;
}


// function to compute the Roes average flux for ith internal face
// takes in comnined vector of left and right conservative variables to compute the flux at Guass point
//                       rhoL rhoUL rhoVL rhoEL rhoR  rhoUR	 rhoVR 	rhoER
//format for arg vector (vl1  vl2   vl3   vl4   vr1   vr2    vr3    vr4) //conservative variable vector
//											 0		1	 	  2		  3	    4		  5		   6		  7	
// takes in states of both left and right and the normal vector to the face its computing the flux at 
matrix2d Flux::Roe::diss2D(std::vector<double> &LR,double &nx, double &ny)
{
	matrix2d dissflux(4,1); //init empty container to store in flux at interface
	//get left and right states
	std::vector<double> Left_state(LR.begin(),LR.begin()+4);
	std::vector<double> Right_state(LR.begin()+4, LR.end());
	//get primitive variables 
	std::vector<double> PrimL = Flux::Roe::Cons2Prim(Left_state); //get primary variables on the left  
	std::vector<double> PrimR = Flux::Roe::Cons2Prim(Right_state); // get primary variables on the right
	
	//get Roes averaged vars;
	std::vector<double> Ravg = Flux::Roe::Roe_avg(PrimL,PrimR);
		// get lamda values
	std::vector<double> lamda = Flux::Roe::wave_speed(Ravg,nx,ny);
	//compute the flux 
	dissflux(0,0);
	dissflux(1,0);
	dissflux(2,0);
	dissflux(3,0);
	return dissflux;

}



	
