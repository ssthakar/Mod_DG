/*
 * 2D DG(Pn) solver using mMatrix2 implementation to store data S S Thakar 27th March 2023
 * */
#ifndef MESH_H
#define MESH_H
#include <chrono>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <numeric>
#include <assert.h>
#include <functional>
#include <algorithm>
#include <cmath>
#include <typeinfo>
#include "mMatrix.h"
#include "mMatrix3.h"


// values that are constant throughout the simulation go here
namespace const_properties
{
	const double gamma = 1.4;
	const double lim_zero = 1e-15;
}
namespace EOS
{
	//equation of state for a monoatomic_gas;
	double perf_gas(double rho, double E,double u,double v); // compute the pressure from the equation of state
}

// convenient typedefs for ease of use
// 2d arrays 
typedef mMatrix2<int> matrix2i;
typedef mMatrix2<double> matrix2d;
typedef mMatrix2<float> matrix2f;
// 3d arrays
typedef mMatrix3<int> matrix3i;
typedef mMatrix3<double> matrix3d;
class reader 
{
	public:
		std::vector<std::vector<double>> readv(std::string s1);
		matrix2d read(std::string s1,int ncols); //read in from .dat file and store in buffer(works for small meshes, needs overhaul for bigger ones)
};

class mesh 
{
	public:
		
		//mesh(int n,int m, int val):mati(m,n,val){};
		// ndimn is to see the dimension of the problem, ntype tells us the kind of cells used in the mesh, nbface tell us the number of boundary faces
		int nelem,ndimn,ntype,npoin,nbface;
		int neqns;
		int ndegr; //taken from control file,// 1 for P(0) 3 for P(1) and 6 for P(2)
		int ngauss_domn; // the number of gasss points for the domain integral taken from control file 
		int ngauss_boun; // the number of gauss points for the boundary of each cell taken from control file
		double A_min; //store in the smallest cell size for CFL condition and timestep calculation
		
		//control will store in inlet velocity, order of polynomial etc
		matrix2d control; //matrix to store in data from control file

		//mesh members
		matrix2i inpoel; // connectivity matrix for mesh 
		matrix2d coords;// coords matrix 
		matrix2i bface; //matrix to store in boundary face and the flags needed to identify a face as inlet,outlet or wall
		matrix2i esup1; //elements surrounding points data struct
		matrix2i esup2; // support data struct for abve
		matrix2i esuel; //element surrounding element data structure
		matrix2i intface; //interface connectivity matrix 
		matrix2d geoface; // data structure to store in normal vector for each face, and gauss points 
		matrix2d geoel;// data structure to store in domain data for each element
		
		//solution containers
		matrix3d U; //3d array that stores in 
		
		matrix3d rhsel; //3d array to store the rhs of every element to solve the time marching 
		
		//for finite volume method the basis function is 1
		// fv_U(i,) = rho | U | V | E 
		matrix2d fv_U; //2d array to store in solution for P(0) DG finite volume

		//mesh methods
		mesh(std::string s1, std::string s2); // constructor reads in mesh file, and control file
		void pre_proc(); //computes all data structures needed in the simulation
		


};

//DG methods 
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


#endif 


