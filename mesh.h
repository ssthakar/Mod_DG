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
	inline double perf_gas(double rho, double E,double u,double v); // compute the pressure from the equation of state
}

// convenient typedefs for ease of use
// 2d arrays 
typedef mMatrix2<int> matrix2i;
typedef mMatrix2<double> matrix2d;
typedef mMatrix2<float> matrix2f;
// 3d arrays
typedef mMatrix3<int> matrix3i;
typedef mMatrix3<double> matrix3d;

namespace reader 
{

	std::vector<std::vector<double>> readv(std::string s1);
	matrix2d read(std::string s1,int ncols); //read in from .dat file and store in buffer(works for small meshes, needs overhaul for bigger ones)
};



//add weights for gauss points
namespace grid
{
	class mesh 
	{
		public:
			int nelem; //number of grid cells without periodic boundary conditions
			int ndimn; // dimension of problem 2D or 3d
			int ntype; // type of cell (does not support hybrid cells) 3 for triangle etc
			int npoin; // number of nodes
			int nbface;// number of boundary faces
			int neqns; // number of conservative variables
			int ndegr; //taken from control file,// 1 for P(0) 3 for P(1) and 6 for P(2)
			int ngauss_domn; // the number of gasss points for the domain integral taken from control file 
			int ngauss_boun; // the number of gauss points for the boundary of each cell taken from control file
			double A_min; //store in the smallest cell size for CFL condition and timestep calculation
			int nmaxface; //total number of faces in the mesh including internal and boundary
			int nintface; // number of internal faces
			//control format ndegr|ngauss_boun|ngauss_domn
			matrix2d control; //matrix to store in data from control file

			//mesh members
			matrix2i inpoel; // connectivity matrix for mesh 
			matrix2d coords;// coords matrix 
			matrix2i bface; //matrix to store in boundary face and the flags needed to identify a face as inlet,outlet or wall
			matrix2i esup1; //elements surrounding points data struct
			matrix2i esup2; // support data struct for abve
			matrix2i esuel; //element surrounding element data structure
			matrix2i intface; //interface connectivity matrix 
			matrix2d geoface; // data structure to store in boundary face data
			matrix2d int_geoface; // data structure to store in internal face data

			//store in jacobian, shape function integrals guass point locations 
			matrix2d domn_weights; // vectors to store in domain weights depending on number of gauss points used 
			matrix2d line_weights; //vector to store in weights for the boundary integral depending on number of guass points used
			matrix2d geoel;// data structure to store in domain data for each element
				
			//solution containers
			matrix3d U; //3d array that stores in the solution unknowns. for each variable
			matrix3d rhsel; //3d array to store rhs for each element
			
			// fv_U(i,) = rho | U | V | E 
			matrix2d fv_U; //2d array to store in solution for P(0) DG finite volume
			
			//mesh methods
			mesh(std::string s1, std::string s2); // constructor reads in mesh file, and control file

	};

	// subroutines to generate data structures for mesh 
	namespace pre_proc
	{
		void set_esup(grid::mesh &mesh1); //generates esup1 and esup2 
		void set_esuel(grid::mesh &mesh1); //generates elsuel
		void set_intfafce(grid::mesh &mesh1); //generate inter face connectivity
		void set_bface(grid::mesh &mesh1);
		void set_boun_geoface(grid::mesh &mesh1); //generates face data 
		void set_geoel(mesh &mesh1); //generates element data
		double el_jacobian(double &x1,double &x2,double &x3,double &y1,double &y2,double &y3); // calculate element jacobian 
		double len(double &x1, double &x2, double &y1, double &y2); // calculate the length of a face 
	}

	namespace post_proc
	{
		void writevtk_mesh(grid::mesh &mesh1,std::string file_name);
	}
}


#endif 


